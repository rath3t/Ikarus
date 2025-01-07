// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "checkfebyautodiff.hh"
#include "resultcollection.hh"
#include "testcommon.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <Eigen/Core>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mechanics/nonlinearelastic.hh>
#include <ikarus/io/resultevaluators.hh>
#include <ikarus/io/resultfunction.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverfactory.hh>
#include <ikarus/solver/nonlinearsolver/trustregion.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/nonlinearoperator.hh>
#include <ikarus/utils/nonlinopfactory.hh>
#include <ikarus/utils/observer/controlvtkwriter.hh>

using Dune::TestSuite;

struct OwnResultFunction
{
  double operator()(const auto& resultArray, [[maybe_unused]] int /* comp */) const { return 7; }
  static std::string name() { return "Seven"; }
  static int ncomps() { return 1; }
};

template <typename Grid, typename Material>
auto NonLinearElasticityLoadControlNRandTR(const Material& mat) {
  TestSuite t("NonLinearElasticityLoadControlNRandTR " + Dune::className(Grid{}) + " " + Dune::className(mat));
  auto grid     = createGrid<Grid>();
  auto gridView = grid->leafGridView();

  using GridView = decltype(gridView);
  using namespace Ikarus;
  using namespace Dune::Functions::BasisFactory;

  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));
  auto vL    = []([[maybe_unused]] auto& globalCoord, auto& lamb) {
    Eigen::Vector2d fext;
    fext.setZero();
    fext[1] = 2 * lamb;
    fext[0] = lamb;
    return fext;
  };

  auto reducedMat = planeStress(mat, 1e-8);
  auto sk         = skills(nonLinearElastic(reducedMat), volumeLoad<2>(vL));
  using FEType    = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;

  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  auto basisP = std::make_shared<const decltype(basis)>(basis);
  Ikarus::DirichletValues dirichletValues(basisP->flat());

  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[1]) < 1e-8)
      dirichletFlags[localView.index(localIndex)] = true;
  });

  auto sparseAssembler = makeSparseFlatAssembler(fes, dirichletValues);

  Eigen::VectorXd d;
  d.setZero(basis.flat().size());
  double lambda = 0.0;

  auto req = typename FEType::Requirement(d, lambda);

  sparseAssembler->bind(req, Ikarus::AffordanceCollections::elastoStatics);
  auto nonLinOp = Ikarus::NonLinearOperatorFactory::op(sparseAssembler, DBCOption::Reduced);

  const double gradTol = 1e-8;

  TrustRegionConfig<> trConfig{
      .parameters = {.verbosity = 1,
                     .maxIter   = 1000,
                     .grad_tol  = gradTol,
                     .corr_tol  = 1e-16, // everything should converge to the gradient tolerance
                     .useRand   = false,
                     .rho_reg   = 1e8,
                     .Delta0    = 1}
  };

  Ikarus::NonlinearSolverFactory trFactory(trConfig);
  auto tr = trFactory.create(sparseAssembler);

  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>>>(
      basis.flat(), d, 2);
  vtkWriter->setFileNamePrefix("Test2DSolid");
  vtkWriter->setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 2);

  auto lc = Ikarus::LoadControl(tr, 1, {0, 50});
  lc.subscribeAll(vtkWriter);
  const auto controlInfo = lc.run();
  nonLinOp.template update<0>();
  const auto maxDisp = std::ranges::max(d);
  double energyExpected;
  if (std::is_same_v<Grid, Grids::Yasp>)
    energyExpected = -2.9605187645668578078;
  else if (std::is_same_v<Grid, Grids::Alu>)
    energyExpected = -2.9530594665063669702;
  else /* std::is_same_v<Grid, Grids::Iga> */
    energyExpected = -1.4533281398929942529;

  double maxDispExpected;
  if (std::is_same_v<Grid, Grids::Yasp>)
    maxDispExpected = 0.11293260007792008115;
  else if (std::is_same_v<Grid, Grids::Alu>)
    maxDispExpected = 0.1123397197762363714;
  else /* std::is_same_v<Grid, Grids::Iga> */
    maxDispExpected = 0.061647849558021668159;

  std::cout << std::setprecision(20) << nonLinOp.value() << std::endl;
  std::cout << "Maxdisp: " << maxDisp << std::endl;
  if constexpr (std::is_same_v<Material, Ikarus::StVenantKirchhoff>) {
    t.check(Dune::FloatCmp::eq(energyExpected, nonLinOp.value()), "energyExpected == nonLinOp.value()")
        << "energyExpected: " << energyExpected << "\nnonLinOp.value(): " << nonLinOp.value();

    t.check(std::abs(maxDispExpected - maxDisp) < 1e-12, "maxDispExpected-maxDisp")
        << "\nmaxDispExpected: \n"
        << maxDispExpected << "\nmaxDisp: \n"
        << maxDisp;
  } else { // using a Neohooke material yields a lower energy and larger displacements
    t.check(Dune::FloatCmp::gt(energyExpected, nonLinOp.value()), "energyExpected > nonLinOp.value()")
        << "energyExpected: " << energyExpected << "\nnonLinOp.value(): " << nonLinOp.value();

    t.check(maxDispExpected < maxDisp, "maxDispExpected<maxDisp") << "maxDispExpected: \n"
                                                                  << maxDispExpected << "\nmaxDisp: \n"
                                                                  << maxDisp;
  }

  Dune::Vtk::VtkWriter<GridView> vtkWriter2(gridView);
  auto resultFunction = makeResultFunction<ResultTypes::PK2Stress>(sparseAssembler);

  t.check(resultFunction->name() == "PK2Stress")
      << "Test resultName: " << resultFunction->name() << "should be PK2Stress";
  t.check(resultFunction->ncomps() == 3) << "Test result comps: " << resultFunction->ncomps() << "should be 3";

  vtkWriter2.addPointData(Dune::Vtk::Function<GridView>(resultFunction));

  auto resultFunction2 =
      makeResultFunction<ResultTypes::PK2Stress, ResultEvaluators::PrincipalStress<2>>(sparseAssembler);
  t.check(resultFunction2->name() == "PrincipalStress")
      << "Test resultName: " << resultFunction2->name() << "should be PrincipalStress";

  t.check(resultFunction2->ncomps() == 2) << "Test result comps: " << resultFunction2->ncomps() << "should be 2";
  vtkWriter2.addPointData(Dune::Vtk::Function<GridView>(resultFunction2));

  auto resultFunction3 = makeResultFunction<ResultTypes::PK2Stress, ResultEvaluators::VonMises>(sparseAssembler);
  t.check(resultFunction3->name() == "VonMises")
      << "Test resultName: " << resultFunction2->name() << "should be VonMises";
  t.check(resultFunction3->ncomps() == 1) << "Test result comps: " << resultFunction2->ncomps() << "should be 1";
  vtkWriter2.addPointData(Dune::Vtk::Function<GridView>(resultFunction3));

  auto resultFunction4 = makeResultVtkFunction<ResultTypes::PK2Stress, OwnResultFunction>(sparseAssembler);
  vtkWriter2.addPointData(resultFunction4);
  vtkWriter2.write("EndResult" + Dune::className<Grid>());

  nonLinOp.template update<1>();
  t.check(controlInfo.success, "Successful result");
  t.check(gradTol >= nonLinOp.derivative().norm(), "Gradient Tolerance should be larger than actual tolerance");
  return t;
}

/**
 * For zero displacements, the stiffness matrix from geometrically linear and non-linear case should be the same,
 * there by ensuring correctness of the Green-Lagrange strain
 */
template <int gridDim, typename Material>
auto GreenLagrangeStrainTest(const Material& mat) {
  static_assert(gridDim == 2 or gridDim == 3, "Only supported for the 2D and 3D case");
  TestSuite t("GreenLagrangeStrainTest with gridDim = " + std::to_string(gridDim));

  const double tol = 1e-8;
  auto grid        = createUGGridFromCorners<gridDim>(CornerDistortionFlag::randomlyDistorted);
  auto gridView    = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis   = Ikarus::makeBasis(gridView, power<gridDim>(lagrange<1>()));
  auto element = gridView.template begin<0>();
  auto nDOF    = basis.flat().size();

  auto fe     = makeFE(basis, skills(nonLinearElastic(mat)));
  auto linMat = [&]() {
    auto linMat = Ikarus::LinearElasticity{mat.materialParameters()};
    if constexpr (gridDim == 3)
      return linMat;
    else {
      if constexpr (isPlaneStress<Material>)
        return Ikarus::planeStress(linMat);
      else
        return Ikarus::planeStrain(linMat);
    }
  }();
  auto feLE = makeFE(basis, skills(Ikarus::linearElastic(linMat)));
  fe.bind(*element);
  feLE.bind(*element);
  Eigen::VectorXd d;
  d.setZero(nDOF);
  double lambda = 0.0;

  auto req   = typename decltype(fe)::Requirement(d, lambda);
  auto reqLE = typename decltype(feLE)::Requirement(d, lambda);

  Eigen::MatrixXd K, KLE;
  K.setZero(nDOF, nDOF);
  KLE.setZero(nDOF, nDOF);

  calculateMatrix(fe, req, Ikarus::MatrixAffordance::stiffness, K);
  calculateMatrix(feLE, reqLE, Ikarus::MatrixAffordance::stiffness, KLE);

  t.check(isApproxSame(K, KLE, tol),
          "Mismatch between linear and non-linear stiffness matrix for zero displacements with gridDim = " +
              std::to_string(gridDim))
      << "\n"
      << K << "\nand \n"
      << KLE << "\n Diff:\n"
      << K - KLE;

  return t;
}

template <typename Material, int gridDim = 2>
auto SingleElementTest(const Material& mat) {
  static_assert(gridDim == 2, "Single element test is applicable only for the 2D case");
  TestSuite t("Single element test for non-linear Q1 element");
  using namespace Ikarus;

  auto grid     = createUGGridFromCorners<gridDim>(CornerDistortionFlag::fixedDistorted);
  auto gridView = grid->leafGridView();

  using namespace Ikarus;

  using namespace Dune::Functions::BasisFactory;
  auto basis       = Ikarus::makeBasis(gridView, power<gridDim>(lagrange<1>()));
  auto element     = gridView.template begin<0>();
  auto nDOF        = basis.flat().size();
  const double tol = 1e-10;

  auto fe = makeFE(basis, skills(nonLinearElastic(mat)));
  fe.bind(*element);

  Eigen::VectorXd d;
  d.setZero(nDOF);
  double lambda = 0.0;

  d << 2, 4, 3.25, -1.2, 0.003, 6, 3, 2.864;

  auto req = typename decltype(fe)::Requirement(d, lambda);

  Eigen::MatrixXd K;
  K.setZero(nDOF, nDOF);
  calculateMatrix(fe, req, Ikarus::MatrixAffordance::stiffness, K);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> essaK;
  essaK.compute(K);
  auto eigenValuesComputed = essaK.eigenvalues();

  /// The eigen values are applicable only if 2x2 Gauss integration points are used
  Eigen::VectorXd eigenValuesExpected;
  eigenValuesExpected.setZero(basis.flat().size());
  eigenValuesExpected << 1e-16, 1e-16, 1845.6296388251504753, 14192.4707553121224317, 19964.32719133414782,
      29973.7943273325380486, 46641.183728849332812, 95447.6156712376251918;
  for (size_t i = 0; i < basis.flat().size(); ++i) {
    if (abs(eigenValuesComputed[i]) > tol) {
      t.check(Dune::FloatCmp::eq(abs(eigenValuesComputed[i]), eigenValuesExpected[i], tol),
              "Mismatch in the " + std::to_string(i + 1) +
                  "-th eigen value in single element test for four node non-linear 2D element");
    }
  }
  return t;
}

template <int gridDim, typename MAT, typename TestSuiteType>
void autoDiffTest(TestSuiteType& t, const MAT& mat, const std::string& testName = "", double tol = 1e-10) {
  using namespace Ikarus;
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

  auto grid     = createUGGridFromCorners<gridDim>(CornerDistortionFlag::randomlyDistorted);
  auto gridView = grid->leafGridView();

  /// We artificially apply a Neumann load on the complete boundary
  Dune::BitSetVector<1> neumannVertices(gridView.size(gridDim), true);
  BoundaryPatch neumannBoundary(gridView, neumannVertices);

  t.subTest(checkFESByAutoDiff(
      gridView, power<gridDim>(lagrange<1>()),
      skills(Ikarus::nonLinearElastic(mat), volumeLoad<gridDim>(vL), neumannBoundaryLoad(&neumannBoundary, nBL)),
      Ikarus::AffordanceCollections::elastoStatics, testName, tol));
}