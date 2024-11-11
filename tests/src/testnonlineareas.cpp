// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
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

#include <Eigen/Core>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/finiteelements/autodifffe.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/nonlinearelastic.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverfactory.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/functionhelper.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/nonlinearoperator.hh>
#include <ikarus/utils/nonlinopfactory.hh>
#include <ikarus/utils/observer/controllogger.hh>
#include <ikarus/utils/observer/controlvtkwriter.hh>
#include <ikarus/utils/observer/nonlinearsolverlogger.hh>

using namespace Ikarus;
using Dune::TestSuite;

template <typename MAT>
auto cantileverBeamTest(const MAT& reducedMat) {
  static_assert(MAT::isReduced, "cantileverBeamTest is only valid for a reduced material (planeStress or planeStrain)");
  TestSuite t("Cantilever Beam for Nonlinear EAS element (Q1E4)");
  constexpr int gridDim = 2;
  using Grid            = Dune::YaspGrid<gridDim>;
  const double L        = 10;
  const double h        = 2;

  Dune::FieldVector<double, gridDim> bbox       = {L, h};
  std::array<int, gridDim> elementsPerDirection = {10, 1};
  auto grid                                     = std::make_shared<Grid>(bbox, elementsPerDirection);
  auto gridView                                 = grid->leafGridView();
  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));

  auto sk      = skills(nonLinearElastic(reducedMat), eas<StrainTags::greenLagrangian>(4));
  using FEType = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;

  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  auto basisP = std::make_shared<const decltype(basis)>(basis);
  Ikarus::DirichletValues dirichletValues(basisP->flat());

  // fix left edge (x=0)
  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[0]) < 1e-8)
      dirichletFlags[localView.index(localIndex)] = true;
  });

  using SparseAssembler = SparseFlatAssembler<decltype(fes), decltype(dirichletValues)>;
  SparseAssembler sparseFlatAssembler(fes, dirichletValues);
  auto sparseAssemblerAM = makeAssemblerManipulator(sparseFlatAssembler);

  Eigen::VectorXd d;
  d.setZero(basis.flat().size());
  double lambda = 0.0;

  auto req = typename FEType::Requirement(d, lambda);

  sparseAssemblerAM->bind(req, Ikarus::AffordanceCollections::elastoStatics, DBCOption::Full);

  // Apply constant point load at the top right corner
  Dune::FieldVector<double, gridDim> topRightPos{L, h};
  const auto globalIndices = utils::globalIndexFromGlobalPosition(basis.flat(), topRightPos);
  auto pointLoad           = [&](const auto&, const auto& par, auto, auto, Eigen::VectorXd& vec) -> void {
    auto loadFactor = par.parameter();
    vec[globalIndices[1]] -= -loadFactor * 1.0;
  };
  sparseAssemblerAM->bind(pointLoad);

  auto linSolver = LinearSolver(SolverTypeTag::sd_UmfPackLU);

  AffordanceCollection elastoStaticsNoScalar(VectorAffordance::forces, MatrixAffordance::stiffness);

  auto nonOp = NonLinearOperatorFactory::op(sparseAssemblerAM, elastoStaticsNoScalar, sparseAssemblerAM->dBCOption());

  constexpr double tol = 1e-10;

  auto nrConfig =
      Ikarus::NewtonRaphsonConfig<decltype(linSolver)>{.parameters = {.tol = tol}, .linearSolver = linSolver};
  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
  auto pathFollowingObserver   = std::make_shared<ControlLogger>();
  auto vtkWriter =
      std::make_shared<ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>>>(basis.flat(), d);
  vtkWriter->setFileNamePrefix("CantileverNonlinearEAS");
  vtkWriter->setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 2);
  auto nr = createNonlinearSolver(nrConfig, nonOp);
  auto lc = LoadControl(nr, 20, {0, 1});
  nr->subscribeAll(nonLinearSolverObserver);
  lc.subscribeAll({pathFollowingObserver, vtkWriter});

  const auto controlInfo = lc.run();

  double expectedLambda  = 1.0;
  double expectedMaxDisp = std::is_same_v<typename MAT::Underlying, NeoHookeT<double>>
                               ? 4.492526443429457
                               : 4.459851990227056; // abs(maxDisp) in ANSYS APDL = 4.48777 (SVK) and = 4.50145 (NH)

  t.check(controlInfo.success);
  const auto maxDisp = std::ranges::max(d.cwiseAbs());

  checkScalars(t, maxDisp, expectedMaxDisp, " Max. displacement", tol);
  checkScalars(t, lambda, expectedLambda, " Load factor", tol);

  return t;
}

template <int gridDim, typename TestSuitType, typename MAT>
void easAutoDiffTest(TestSuitType& t, const MAT& mat) {
  using namespace Dune::Functions::BasisFactory;
  std::array<int, (gridDim == 2) ? 4 : 3> easParameters;
  if constexpr (gridDim == 2)
    easParameters = {0, 4, 5, 7};
  else
    easParameters = {0, 9, 21};

  auto grid = createUGGridFromCorners<gridDim>(CornerDistortionFlag::randomlyDistorted);
  grid->globalRefine(2);
  auto gridView = grid->leafGridView();

  for (const int numberOfEASParameters : easParameters) {
    auto sk = skills(nonLinearElastic(mat), eas<StrainTags::greenLagrangian>(numberOfEASParameters));
    t.subTest(checkFESByAutoDiff(gridView, power<gridDim>(lagrange<1>()), sk, AffordanceCollections::elastoStatics,
                                 " (numberOfEASParameters = " + std::to_string(numberOfEASParameters) + ")"));
  }
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t("Nonlinear EAS Test");
  auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = 100.0, .nu = 0.3});
  StVenantKirchhoff matSVK(matParameter);
  NeoHooke matNH(matParameter);
  auto reducedMatSVK = planeStrain(matSVK);
  auto reducedMatNH  = planeStrain(matNH);

  easAutoDiffTest<2>(t, reducedMatSVK);
  easAutoDiffTest<3>(t, matSVK);

  easAutoDiffTest<2>(t, reducedMatNH);
  easAutoDiffTest<3>(t, matNH);

  t.subTest(cantileverBeamTest(reducedMatSVK));
  t.subTest(cantileverBeamTest(reducedMatNH));
  return t.exit();
}