// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "checkfebyautodiff.hh"
#include "testcommon.hh"
#include "testhelpers.hh"

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/truss.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverfactory.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controlvtkwriter.hh>
#include <ikarus/utils/observer/genericobserver.hh>
#include <ikarus/utils/observer/nonlinearsolverlogger.hh>

using namespace Ikarus;
using Dune::TestSuite;

Eigen::Matrix<double, 4, 4> rotationMatrix2D(double theta) {
  const auto c = std::cos(theta);
  const auto s = std::sin(theta);
  Eigen::Matrix2d rotMat;
  rotMat.setZero();
  rotMat << c, s, -s, c;
  Eigen::Matrix<double, 4, 4> Q;
  Q.setZero();
  Q.block<2, 2>(0, 0) = rotMat;
  Q.block<2, 2>(2, 2) = rotMat;
  return Q;
}

static auto vonMisesTrussTest() {
  TestSuite t("vonMisesTrussTest");
  /// Construct grid
  constexpr double h = 1.0;
  constexpr double L = 10.0;
  auto grid          = createGrid<Grids::OneDFoamGridIn2D>();
  auto gridView      = grid->leafGridView();

  /// Construct basis
  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>()));

  /// Create finite elements
  constexpr double E   = 30000;
  constexpr double A   = 0.1;
  constexpr double tol = 1e-10;
  auto sk              = skills(truss(E, A));
  using FEType         = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;
  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  /// Collect dirichlet nodes
  auto basisP = std::make_shared<const decltype(basis)>(basis);
  DirichletValues dirichletValues(basisP->flat());
  dirichletValues.fixBoundaryDOFs(
      [&](auto& dirichletFlags, auto&& globalIndex) { dirichletFlags[globalIndex] = true; });
  t.check(dirichletValues.fixedDOFsize() == 4, "Number of fixed DOFs is not equal to 4");

  /// Create assembler
  auto denseFlatAssembler = makeAssemblerManipulator(DenseFlatAssembler(fes, dirichletValues));

  /// Create non-linear operator
  double lambda = 0.0;
  Eigen::VectorXd d;
  d.setZero(basis.flat().size());

  auto req = FEType::Requirement();
  req.insertGlobalSolution(d).insertParameter(lambda);
  denseFlatAssembler->bind(req, AffordanceCollections::elastoStatics, DBCOption::Full);

  auto pointLoad = [&](const auto&, const auto&, auto, auto, Eigen::VectorXd& vec) -> void { vec[3] -= -lambda; };
  denseFlatAssembler->bind(pointLoad);

  /// Test tangent stiffness matrix for element 0 for geometrically linear case
  const auto& ele0 = fes[0].gridElement();
  Eigen::MatrixXd ke0;
  ke0.setZero(4, 4);
  calculateMatrix(fes[0], req, MatrixAffordance::stiffness, ke0);
  const auto X1            = Dune::toEigen(ele0.geometry().corner(0));
  const auto X2            = Dune::toEigen(ele0.geometry().corner(1));
  const Eigen::Vector2d A1 = (X2 - X1).eval();
  const Eigen::Vector2d A0 = {1, 0}; // horizontal axis
  const double L0          = sqrt(A1.squaredNorm());
  const double theta0      = std::acos(A1.dot(A0) / (A0.norm() * A1.norm()));
  Eigen::Matrix<double, 4, 4> kbar;
  kbar.setZero();
  kbar(0, 0) = kbar(2, 2) = 1.0;
  kbar(0, 2) = kbar(2, 0) = -1.0;
  const auto Q            = rotationMatrix2D(theta0);
  const auto k0           = (Q.transpose() * kbar * Q) * E * A / L0;
  t.check(isApproxSame(k0, ke0, tol), "(k0 is not equal to ke0 for geometrically linear case.") << "Expected:\n"
                                                                                                << k0 << "\nActual:\n"
                                                                                                << ke0;

  /// Choose linear solver
  auto linSolver = LinearSolver(SolverTypeTag::d_LDLT);

  NewtonRaphsonConfig<decltype(linSolver)> nrConfig{.linearSolver = linSolver};

  NonlinearSolverFactory nrFactory(nrConfig);
  auto nr = nrFactory.create(denseFlatAssembler);

  /// Create Observer to write information of the non-linear solver
  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
  auto nonLinOp                = Ikarus::NonLinearOperatorFactory::op(denseFlatAssembler);

  t.check(utils::checkGradient(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << "Check gradient failed";
  t.check(utils::checkHessian(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << "Check Hessian failed";

  constexpr int loadSteps = 10;

  Eigen::Matrix3Xd lambdaAndDisp;
  lambdaAndDisp.setZero(Eigen::NoChange, loadSteps + 1);
  /// Create Observer which executes when control routines messages
  auto lvkObserver =
      std::make_shared<GenericObserver<ControlMessages>>(ControlMessages::SOLUTION_CHANGED, [&](int step) {
        lambdaAndDisp(0, step) = lambda; // load factor
        lambdaAndDisp(1, step) = d[2];   // horizontal displacement at center node
        lambdaAndDisp(2, step) = d[3];   // vertical displacement at center node
      });

  /// Create Observer which writes vtk files when control routines messages
  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>>>(
      basis.flat(), d, 2);
  vtkWriter->setFieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2);
  vtkWriter->setFileNamePrefix("vonMisesTruss");

  /// Create loadcontrol
  auto lc = LoadControl(nr, loadSteps, {0, 0.5});
  lc.nonlinearSolver().subscribeAll(nonLinearSolverObserver);
  lc.subscribeAll({vtkWriter, lvkObserver});

  /// Execute!
  auto controlInfo = lc.run();
  t.check(controlInfo.success == true) << "Load control failed to converge";

  Eigen::VectorXd lambdaVec = lambdaAndDisp.row(0);
  Eigen::VectorXd uVec      = lambdaAndDisp.row(1);
  Eigen::VectorXd vVec      = -lambdaAndDisp.row(2);

  // return the load factor as a function of the vertical displacement
  auto analyticalLoadDisplacementCurve = [&](auto& v) {
    const double Ltruss = std::sqrt(h * h + L * L);
    return 2.0 * E * A * Dune::power(h, 3) / Dune::power(Ltruss, 3) *
           (v / h - 1.5 * Dune::power(v / h, 2) + 0.5 * Dune::power(v / h, 3));
  };

  for (std::size_t i = 0; i < lambdaAndDisp.cols(); ++i) {
    checkScalars(t, lambdaVec[i], analyticalLoadDisplacementCurve(vVec[i]), " Incorrect load factor", tol);
    checkScalars(t, uVec[i], 0.0, " Incorrect horizontal displacement", tol);
  }

  double expectedCauchyAxialForce = -2.785092479363964262;
  double expectedPK2AxialForce    = -2.787684078189828;

  // due to the symmetry of the problem, same axial forces are expected in both elements
  for (const auto& fe : fes) {
    const auto N    = fe.calculateAt<ResultTypes::cauchyAxialForce>(req, {0.5}).asVec()[0];
    const auto NPK2 = fe.calculateAt<ResultTypes::PK2AxialForce>(req, {0.5}).asVec()[0];
    checkScalars(t, N, expectedCauchyAxialForce, " Incorrect Cauchy Axial force", tol);
    checkScalars(t, NPK2, expectedPK2AxialForce, " Incorrect PK2 Axial force", tol);
  }

  return t;
}

static auto truss3dTest() {
  TestSuite t("Truss3DTest");
  /// Construct grid
  auto grid     = createGrid<Grids::OneDFoamGridIn3D>();
  auto gridView = grid->leafGridView();

  /// Construct basis
  using namespace Dune::Functions::BasisFactory;
  auto basis      = Ikarus::makeBasis(gridView, power<3>(lagrange<1>()));
  using LocalView = std::remove_cvref_t<decltype(basis.flat().localView())>;

  /// Create finite elements
  constexpr double E   = 2000;
  constexpr double A   = 1000;
  constexpr double tol = 1e-10;
  auto sk              = skills(truss(E, A));
  using FEType         = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;
  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  /// Collect dirichlet nodes
  auto basisP = std::make_shared<const decltype(basis)>(basis);
  DirichletValues dirichletValues(basisP->flat());
  dirichletValues.fixBoundaryDOFs(
      [&](auto& dirichletFlags, auto&& globalIndex) { dirichletFlags[globalIndex] = true; });
  std::array<LocalView::MultiIndex, 1> fixDOFs = {3}; // array containing the DOFs to be fixed at the center node
  for (const auto dof : fixDOFs)
    dirichletValues.setSingleDOF(dof, true);
  t.check(dirichletValues.fixedDOFsize() == 7, "Number of fixed DOFs is not equal to 7");

  /// Create assembler
  auto denseFlatAssembler = makeDenseFlatAssembler(fes, dirichletValues);

  double lambda = 1.0;
  Eigen::VectorXd d;
  d.setZero(basis.flat().size());

  auto req = FEType::Requirement();
  req.insertGlobalSolution(d).insertParameter(lambda);
  denseFlatAssembler->bind(req, AffordanceCollections::elastoStatics, DBCOption::Full);
  const auto& K = denseFlatAssembler->matrix();
  auto R        = denseFlatAssembler->vector();
  R[4]          = lambda;

  /// Choose linear solver
  auto linSolver = LinearSolver(SolverTypeTag::d_LDLT);
  linSolver.compute(K);
  linSolver.solve(d, -R);

  /// Compute eigenvalues of the stiffness matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> essaK;
  essaK.compute(K);
  auto eigenValuesComputed = essaK.eigenvalues();

  Eigen::Vector2d axialForces = Eigen::Vector2d::Zero();
  for (int i = 0; const auto& fe : fes) {
    axialForces[i] = fe.calculateAt<ResultTypes::linearAxialForce>(req, {0.5}).asVec()[0];
    ++i;
  }

  Eigen::VectorXd expectedDisplacements, expectedEigenValues, expectedAxialForces;
  expectedDisplacements.setZero(basis.flat().size());
  expectedEigenValues.setZero(basis.flat().size());
  expectedAxialForces.setZero(2);
  expectedDisplacements << 0, 0, 0, 0, -0.0008521190304919996, 0.0002491175412207767, 0, 0, 0;
  expectedEigenValues << 1, 1, 1, 1, 1, 1, 1, 1081.004736091822, 734025.9250446925;
  expectedAxialForces << 13.78521130992862, -14.91038752769381;

  t.check(isApproxSame(expectedEigenValues, eigenValuesComputed, tol))
      << std::setprecision(16) << "Incorrect eigenvalues of K - Expected: " << expectedEigenValues.transpose()
      << " Actual: " << eigenValuesComputed.transpose();

  t.check(isApproxSame(expectedDisplacements, d, tol))
      << std::setprecision(16) << "Incorrect displacements - Expected: " << expectedDisplacements.transpose()
      << " Actual: " << d.transpose();

  t.check(isApproxSame(expectedAxialForces, axialForces, tol))
      << std::setprecision(16) << "Incorrect axial forces - Expected: " << expectedAxialForces.transpose()
      << " Actual: " << axialForces.transpose();

  return t;
}

template <int worldDim, typename GridType>
auto trussAutoDiffTest() {
  TestSuite t("Truss autodiff");
  using namespace Dune::Functions::BasisFactory;
  {
    auto grid     = createGrid<GridType>();
    auto gridView = grid->leafGridView();
    t.subTest(checkFESByAutoDiff(gridView, power<worldDim>(lagrange<1>()), skills(truss(1000.0, 0.0956)),
                                 AffordanceCollections::elastoStatics));
  }
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t("TrussTest");
  t.subTest(trussAutoDiffTest<2, Grids::OneDFoamGridIn2D>());
  t.subTest(trussAutoDiffTest<3, Grids::OneDFoamGridIn3D>());
  t.subTest(vonMisesTrussTest());
  t.subTest(truss3dTest());
  return t.exit();
}
