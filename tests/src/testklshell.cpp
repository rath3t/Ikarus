// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "checkfebyautodiff.hh"
#include "testcommon.hh"
#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#if HAVE_DUNE_IGA
  #include <dune/iga/nurbsbasis.hh>
#endif
#include "spdlog/spdlog.h"

#include <Eigen/Core>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/finiteelements/mechanics/kirchhoffloveshell.hh>
#include <ikarus/io/resultfunction.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/solver/nonlinearsolver/trustregion.hh>
#include <ikarus/utils/algorithms.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/nonlinearoperator.hh>
#include <ikarus/utils/nonlinopfactory.hh>
#include <ikarus/utils/observer/controlvtkwriter.hh>
#include <ikarus/utils/observer/nonlinearsolverlogger.hh>

using Dune::TestSuite;

static auto NonLinearKLShellLoadControlTR() {
  TestSuite t("NonLinearKLShellLoadControlTR");
  constexpr auto dimworld        = 3;
  const std::array<int, 2> order = {1, 1};

  const std::array<std::vector<double>, 2> knotSpans = {
      {{0, 0, 1, 1}, {0, 0, 1, 1}}
  };

  using ControlPoint = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints = {
      {{.p = {0, 0, 0}, .w = 1}, {.p = {10, 0, 0}, .w = 1}},
      {{.p = {0, 2, 0}, .w = 1}, {.p = {10, 2, 0}, .w = 1}}
  };

  std::array<int, 2> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<2, dimworld>;

  Dune::IGA::NURBSPatchData<2, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  for (int i = 0; i < 2; ++i)
    patchData = degreeElevate(patchData, i, 1);

  auto grid = std::make_shared<Grid>(patchData);
  grid->globalRefine(2);
  auto gridView = grid->leafGridView();

  using namespace Ikarus;
  using namespace Dune::Functions::BasisFactory;
  const double E         = 1000;
  const double nu        = 0.0;
  const double thickness = 0.1;
  auto basis             = Ikarus::makeBasis(gridView, power<3>(nurbs(), FlatInterleaved()));
  auto vL                = [thickness]([[maybe_unused]] auto& globalCoord, auto& lamb) {
    Eigen::Vector3d fext;
    fext.setZero();
    fext[2] = 2 * Dune::power(thickness, 3) * lamb;
    return fext;
  };

  auto sk = skills(kirchhoffLoveShell({.youngs_modulus = E, .nu = nu, .thickness = thickness}), volumeLoad<3>(vL));

  using FEType = decltype(Ikarus::makeFE(basis, sk));
  std::vector<FEType> fes;

  for (auto&& element : elements(gridView)) {
    fes.emplace_back(Ikarus::makeFE(basis, sk));
    fes.back().bind(element);
  }

  auto basisP = std::make_shared<const decltype(basis)>(basis);
  Ikarus::DirichletValues dirichletValues(basisP->flat());

  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[0]) < 1e-8)
      dirichletFlags[localView.index(localIndex)] = true;
  });

  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[0]) > 10.0 - 1e-8)
      dirichletFlags[localView.index(localIndex)] = true;
  });

  auto sparseAssembler = makeSparseFlatAssembler(fes, dirichletValues);

  Eigen::VectorXd d;
  d.setZero(basis.flat().size());
  double lambda = 0.0;

  auto req = FEType::Requirement();
  req.insertGlobalSolution(d).insertParameter(lambda);

  [[maybe_unused]] auto req2 = FEType::Requirement(d, lambda);

  sparseAssembler->bind(req, Ikarus::AffordanceCollections::elastoStatics);
  auto nonLinOp = Ikarus::NonLinearOperatorFactory::op(sparseAssembler);

  t.check(utils::checkGradient(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << "Check gradient failed";
  t.check(utils::checkHessian(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << "Check Hessian failed";

  const double gradTol = 1e-14;

  auto tr = makeTrustRegion(nonLinOp);
  tr->setup({.verbosity = 1,
             .maxIter   = 1000,
             .grad_tol  = gradTol,
             .corr_tol  = 1e-16, // everything should converge to the gradient tolerance
             .useRand   = false,
             .rho_reg   = 1e8,
             .Delta0    = 1});

  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>>>(
      basis.flat(), d, 2);
  vtkWriter->setFileNamePrefix("TestKLShell");
  vtkWriter->setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 3);

  auto lc = LoadControl(tr, 1, {0, 1});
  lc.subscribeAll(vtkWriter);
  const auto controlInfo = lc.run(req);

  t.check(controlInfo.success);

  const auto maxDisp = std::ranges::max(d);
  std::cout << std::setprecision(16) << maxDisp << std::endl;
  t.check(Dune::FloatCmp::eq(0.2087577577946809, maxDisp, 1e-6))
      << std::setprecision(16) << "The maximum displacement is " << maxDisp << "but it should be " << 0.2087577577946809
      << ". The difference is " << 0.2087577577946809 - maxDisp;
  return t;
}

auto singleElementTest() {
  TestSuite t("Kirchhoff-Love autodiff");
  using namespace Dune::Functions::BasisFactory;

  auto vL = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    Eigen::Vector<typename VectorType::field_type, VectorType::dimension> fExt;
    fExt.setZero();
    fExt[1] = 2 * lamb;
    return fExt;
  };

  auto nBL = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    Eigen::Vector<typename VectorType::field_type, VectorType::dimension> fExt;
    fExt.setZero();
    fExt[0] = lamb / 40;
    return fExt;
  };
  {
    constexpr auto dimworld        = 3;
    const std::array<int, 2> order = {1, 1};

    const std::array<std::vector<double>, 2> knotSpans = {
        {{0, 0, 1, 1}, {0, 0, 1, 1}}
    };

    using ControlPoint = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointType;

    const std::vector<std::vector<ControlPoint>> controlPoints = {
        {{.p = {0, 0, 0}, .w = 1}, {.p = {10, 0, 0}, .w = 1}},
        {{.p = {0, 2, 0}, .w = 1}, {.p = {10, 2, 0}, .w = 1}}
    };

    std::array<int, 2> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};

    auto controlNet = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointNetType(dimsize, controlPoints);
    using Grid      = Dune::IGA::NURBSGrid<2, dimworld>;

    Dune::IGA::NURBSPatchData<2, dimworld> patchData;
    patchData.knotSpans     = knotSpans;
    patchData.degree        = order;
    patchData.controlPoints = controlNet;
    for (int i = 0; i < 2; ++i)
      patchData = degreeElevate(patchData, i, 1);

    auto grid     = std::make_shared<Grid>(patchData);
    auto gridView = grid->leafGridView();
    /// We artificially apply a Neumann load on the complete boundary
    Dune::BitSetVector<1> neumannVertices(gridView.size(2), true);
    BoundaryPatch neumannBoundary(gridView, neumannVertices);

    auto klShell = Ikarus::kirchhoffLoveShell({.youngs_modulus = 1000, .nu = 0.0, .thickness = 0.1});

    t.subTest(checkFESByAutoDiff(
        gridView, power<3>(nurbs()),
        Ikarus::skills(klShell, Ikarus::volumeLoad<3>(vL), Ikarus::neumannBoundaryLoad(&neumannBoundary, nBL)),
        Ikarus::AffordanceCollections::elastoStatics));
  }
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t("Kirchhoff-Love");
  t.subTest(singleElementTest());
  t.subTest(NonLinearKLShellLoadControlTR());
  return t.exit();
}
