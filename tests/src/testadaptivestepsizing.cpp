// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "checkfebyautodiff.hh"
#include "testcommon.hh"
#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/iga/nurbsbasis.hh>

#include <Eigen/Core>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/pathfollowing.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/kirchhoffloveshell.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphsonwithscalarsubsidiaryfunction.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/nonlinearoperator.hh>
#include <ikarus/utils/nonlinopfactory.hh>
#include <ikarus/utils/observer/controllogger.hh>
#include <ikarus/utils/observer/controlvtkwriter.hh>
#include <ikarus/utils/observer/nonlinearsolverlogger.hh>

using Dune::TestSuite;

template <typename PathFollowingType>
auto KLShellAndAdaptiveStepSizing(const PathFollowingType& pft, const std::vector<std::vector<int>>& expectedIterations,
                                  const std::vector<std::vector<double>>& expectedResults, const int targetIterations,
                                  const double stepSize) {
  TestSuite t("KLShellAndAdaptiveStepSizing --> " + pft.name());
  constexpr auto dimWorld        = 3;
  const std::array<int, 2> order = {2, 1};

  const std::array<std::vector<double>, 2> knotSpans = {
      {{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}
  };

  using ControlPoint = Dune::IGA::NURBSPatchData<2, dimWorld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints = {
      {   {.p = {-0.5, 0, 0.0}, .w = 10},     {.p = {-0.5, 1.0, 0.0}, .w = 1}},
      {{.p = {0, 0, 0.2}, .w = 0.766044}, {.p = {0, 1.0, 0.2}, .w = 0.766044}},
      {     {.p = {0.5, 0, 0.0}, .w = 1},      {.p = {0.5, 1.0, 0.0}, .w = 1}}
  };

  std::array<int, 2> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<2, dimWorld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<2, dimWorld>;

  Dune::IGA::NURBSPatchData<2, dimWorld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  using namespace Dune::IGA;
  patchData = degreeElevate(patchData, 1, 1);

  auto grid = std::make_shared<Grid>(patchData);

  grid->globalRefineInDirection(0, 2);
  grid->globalRefineInDirection(1, 2);
  auto gridView        = grid->leafGridView();
  const auto& indexSet = gridView.indexSet();

  using GridView = decltype(gridView);
  using namespace Ikarus;
  using namespace Dune::Functions::BasisFactory;
  const double E         = 4.32e4;
  const double nu        = 0.3;
  const double thickness = 0.00001;
  auto basis             = Ikarus::makeBasis(gridView, power<3>(nurbs()));

  /// Neumann boundary condition at right edge (x=0.5)
  Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
  for (auto&& vertex : vertices(gridView))
    if (std::abs(vertex.geometry().corner(0)[0] - 0.5) < 1e-8)
      neumannVertices[indexSet.index(vertex)] = true;

  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);

  auto nBL = [&](auto& globalCoord, auto& lamb) {
    Eigen::Vector3d fext;
    fext.setZero();
    fext[0] = -lamb;
    return fext;
  };

  Ikarus::KlArgs args = {E, nu, thickness};
  t.subTest(checkFESByAutoDiff(gridView, power<3>(nurbs()),
                               skills(kirchhoffLoveShell(args), neumannBoundaryLoad(&neumannBoundary, nBL)),
                               Ikarus::AffordanceCollections::elastoStatics));
  auto sk = skills(kirchhoffLoveShell(args), neumannBoundaryLoad(&neumannBoundary, nBL));

  using FEType = decltype(Ikarus::makeFE(basis, sk));
  std::vector<FEType> fes;

  for (auto& element : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(element);
  }

  auto basisP = std::make_shared<const decltype(basis)>(basis);
  DirichletValues dirichletValues(basisP->flat());

  /// fix all DOFs at the left edge (x=-0.5)
  dirichletValues.fixDOFs([&](auto& basis_, auto&& dirichletFlags) {
    Dune::Functions::forEachBoundaryDOF(basis_, [&](auto&& localIndex, auto&& localView, auto&& intersection) {
      if (std::abs(intersection.geometry().center()[0] + 0.5) < 1e-8)
        dirichletFlags[localView.index(localIndex)] = true;
    });
  });

  /// fix y and z DOFs at the right edge (x=0.5)
  dirichletValues.fixDOFs([&](auto& basis_, auto&& dirichletFlags) {
    for (auto fixedDirection = 1; fixedDirection < 3; ++fixedDirection) {
      Dune::Functions::forEachBoundaryDOF(Dune::Functions::subspaceBasis(basis_, fixedDirection),
                                          [&](auto&& localIndex, auto&& localView, auto&& intersection) {
                                            if (std::abs(intersection.geometry().center()[0] - 0.5) < 1e-8)
                                              dirichletFlags[localView.index(localIndex)] = true;
                                          });
    }
  });

  auto sparseAssembler = makeSparseFlatAssembler(fes, dirichletValues);

  Eigen::VectorXd d;
  d.setZero(basis.flat().size());
  double lambda = 0.0;

  auto req = fes[0].createRequirement();
  req.insertGlobalSolution(d).insertParameter(lambda);
  sparseAssembler->bind(req, Ikarus::AffordanceCollections::elastoStatics);

  auto linSolver = LinearSolver(SolverTypeTag::sd_SimplicialLDLT);

  int loadSteps = 6;

  auto nrConfig = NewtonRaphsonWithSubsidiaryFunctionConfig<decltype(linSolver)>{.linearSolver = linSolver};

  NonlinearSolverFactory nrFactory(nrConfig);
  auto nr  = nrFactory.create(sparseAssembler);
  auto nr2 = nrFactory.create(sparseAssembler);

  auto dass = AdaptiveStepSizing::IterationBased{};
  auto nass = AdaptiveStepSizing::NoOp{};
  dass.setTargetIterations(targetIterations);

  /// control routine with and without step sizing
  auto crWSS  = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft, dass);
  auto crWoSS = Ikarus::PathFollowing(nr2, loadSteps, stepSize, pft, nass);

  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
  auto pathFollowingObserver   = std::make_shared<ControlLogger>();
  crWSS.nonlinearSolver().subscribeAll(nonLinearSolverObserver);
  crWoSS.nonlinearSolver().subscribeAll(nonLinearSolverObserver);

  t.checkThrow<Dune::InvalidStateException>(
      [&]() { nonLinearSolverObserver->update(Ikarus::NonLinearSolverMessages::BEGIN); },
      "nonLinearSolverObserver should have failed for the BEGIN message");

  t.checkThrow<Dune::InvalidStateException>(
      [&]() { nonLinearSolverObserver->update(Ikarus::NonLinearSolverMessages::END); },
      "nonLinearSolverObserver should have failed for the END message");

  /// Create Observer which writes vtk files when control routines messages
  /// SOLUTION_CHANGED
  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>>>(
      basis.flat(), d, 2);
  vtkWriter->setFieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 3);
  vtkWriter->setFileNamePrefix("testAdaptiveStepSizing" + pft.name());

  crWoSS.subscribeAll({vtkWriter, pathFollowingObserver});
  crWoSS.unSubscribeAll(vtkWriter);
  crWSS.subscribeAll({pathFollowingObserver});
  crWSS.subscribe(ControlMessages::SOLUTION_CHANGED, vtkWriter);

  const std::string& message1 = " --> " + pft.name() + " with default adaptive step sizing";
  const std::string& message2 = " --> " + pft.name() + " without default adaptive step sizing";

  t.checkThrow<Dune::InvalidStateException>(
      [&]() {
        auto dass2                 = AdaptiveStepSizing::IterationBased{};
        auto nr3                   = nrFactory.create(sparseAssembler);
        auto crWSS2                = Ikarus::PathFollowing(nr3, loadSteps, stepSize, pft, dass2);
        const auto controlInfoWSS2 = crWSS2.run();
      },
      "IterationBased should fail for targetIterations being 0");

  resetNonLinearOperatorParametersToZero(crWSS.nonlinearSolver().nonLinearOperator());
  const auto controlInfoWSS = crWSS.run();
  const double tolDisp      = 1e-13;
  const double tolLoad      = 1e-12;
  checkScalars(t, std::ranges::max(d), expectedResults[0][0], message1 + " <Max Displacement>", tolDisp);
  checkScalars(t, lambda, expectedResults[0][1], message1 + " <Lambda>", tolLoad);
  resetNonLinearOperatorParametersToZero(crWSS.nonlinearSolver().nonLinearOperator());

  const auto controlInfoWoSS = crWoSS.run();

  checkScalars(t, std::ranges::max(d), expectedResults[1][0], message2 + " <Max Displacement>", tolDisp);
  checkScalars(t, lambda, expectedResults[1][1], message2 + " <Lambda>", tolLoad);
  resetNonLinearOperatorParametersToZero(crWSS.nonlinearSolver().nonLinearOperator());

  const int controlInfoWSSIterations =
      std::accumulate(controlInfoWSS.solverInfos.begin(), controlInfoWSS.solverInfos.end(), 0,
                      [](int a, auto& b) { return b.iterations + a; });

  t.check(controlInfoWSS.success, "No convergence" + message1);
  t.check(controlInfoWoSS.success, "No convergence" + message2);

  t.check(controlInfoWSS.totalIterations < controlInfoWoSS.totalIterations)
      << "Total iterations should be less --> " << controlInfoWSS.totalIterations << " > "
      << std::to_string(controlInfoWoSS.totalIterations) + " --> " + pft.name();

  t.check(controlInfoWSSIterations == controlInfoWSS.totalIterations)
      << "Total number of iterations is wrong --> " << controlInfoWSS.totalIterations << " is not equal to "
      << std::to_string(controlInfoWSSIterations) + " --> " + pft.name();

  checkSolverInfos(t, expectedIterations[0], controlInfoWSS, loadSteps, message1);
  checkSolverInfos(t, expectedIterations[1], controlInfoWoSS, loadSteps, message2);

  auto nonLinOp = NonLinearOperatorFactory::op(sparseAssembler);
  t.check(utils::checkGradient(nonLinOp, {.draw = false})) << "Check gradient failed";
  t.check(utils::checkHessian(nonLinOp, {.draw = false})) << "Check hessian failed";

  return t;
}
#include <cfenv>
int main(int argc, char** argv) {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
  TestSuite t("Adaptive Stepsizing");
  Ikarus::init(argc, argv);

  auto alc = Ikarus::ArcLength{};
  auto lc  = Ikarus::LoadControlSubsidiaryFunction{};

  /// expected iterations for each step for a path following type with and without step sizing
  const std::vector<std::vector<int>> expectedIterationsALC = {
      {9, 9, 6, 5, 4, 4},
      {9, 8, 6, 5, 5, 5}
  };
  const std::vector<std::vector<int>> expectedIterationsLC = {
      {10, 4, 4, 3, 3, 3},
      {10, 6, 5, 4, 4, 4}
  };

  /// expected results(max(displacement) and lambda) for a path following type with and without step sizing
  const std::vector<std::vector<double>> expectedResultsALC = {
      {0.1032139637288574, 0.0003103004514250302},
      { 0.162759603260405, 0.0007765975850229621}
  };
  const std::vector<std::vector<double>> expectedResultsLC = {
      {0.08741028329554552, 0.0002318693543601816},
      {  0.144353999993308,                  6e-4}
  };

  t.subTest(KLShellAndAdaptiveStepSizing(alc, expectedIterationsALC, expectedResultsALC, 3, 0.4));
  t.subTest(KLShellAndAdaptiveStepSizing(lc, expectedIterationsLC, expectedResultsLC, 2, 1e-4));
  return t.exit();
}
