
#include <config.h>
#include <catch2/catch_all.hpp>

#include "testHelpers.hh"
#include "common.hh"

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include "spdlog/spdlog.h"

#include <Eigen/Core>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlRoutines/loadControl.hh>
#include <ikarus/finiteElements/mechanics/nonLinearElasticityFE.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/solver/nonLinearSolver/trustRegion.hh>
#include <ikarus/utils/algorithms.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>

//template <typename T>
//class NonLinearElasticityLoadControlNRandTR : public testing::Test {
//public:
//  using GridId = T;
//  NonLinearElasticityLoadControlNRandTR() : value_{createGrid<T>()} {}
//  decltype(createGrid<T>()) value_;
//};
//using GridTypes = ::testing::Types<Grids::Yasp, Grids::Alu, Grids::Iga>;
//
//TYPED_TEST_SUITE(NonLinearElasticityLoadControlNRandTR, GridTypes);

TEMPLATE_TEST_CASE("NonLinearElasticityLoadControlNRandTR: ComputeMaxDisp", "[nonLinearElasticityTest.cpp]",Grids::Yasp, Grids::Alu, Grids::Iga) {
  auto grid = createGrid<TestType>();
  auto gridView = grid->leafGridView();

  using namespace Ikarus;

  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));

  auto localView = basis.localView();
  std::vector<Ikarus::NonLinearElasticityFE<decltype(basis)>> fes;
  auto volumeLoad = [](auto& globalCoord, auto& lamb) {
    Eigen::Vector2d fext;
    fext.setZero();
    fext[1] = 2 * lamb;
    fext[0] = lamb;
    return fext;
  };
  for (auto& element : elements(gridView))
    fes.emplace_back(basis, element, 1000, 0.3, volumeLoad);

  std::vector<bool> dirichletFlags(basis.size(), false);

  Dune::Functions::forEachBoundaryDOF(basis, [&](auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[1]) < 1e-8) {
      dirichletFlags[localView.index(localIndex)[0]] = true;
    }
  });

  auto sparseAssembler = SparseFlatAssembler(basis, fes, dirichletFlags);

  Eigen::VectorXd d;
  d.setZero(basis.size());
  double lambda = 0.0;

  auto residualFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::VectorAffordances::forces)
                                     .build();
    return sparseAssembler.getVector(req);
  };

  auto KFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::MatrixAffordances::stiffness)
                                     .build();
    return sparseAssembler.getMatrix(req);
  };

  auto energyFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::ScalarAffordances::mechanicalPotentialEnergy)
                                     .build();
    return sparseAssembler.getScalar(req);
  };

  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(energyFunction, residualFunction, KFunction),
                                            parameter(d, lambda));

  const double gradTol = 1e-8;

  auto tr = Ikarus::makeTrustRegion(nonLinOp);
  tr->setup({.verbosity = 1,
             .maxiter   = 1000,
             .grad_tol  = gradTol,
             .corr_tol  = 1e-16,  // everything should converge to the gradient tolerance
             .useRand   = false,
             .rho_reg   = 1e8,
             .Delta0    = 1});

  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<decltype(basis)>>(basis, d, 2);
  vtkWriter->setFileNamePrefix("Test2Dsolid");
  vtkWriter->setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 2);

  auto lc = Ikarus::LoadControl(tr, 1, {0, 2000});
  lc.subscribeAll(vtkWriter);
  const auto controlInfo = lc.run();
  nonLinOp.template update<0>();
  const auto maxDisp = std::ranges::max(d);
  const double energyExpected =  (std::is_same_v<TestType, Grids::Yasp>) ? -1.4809559783564966e+03 : ((std::is_same_v<TestType, Grids::Alu>)  ? -1.4842107484533601e+03 :  /* std::is_same_v<TestType, Grids::Iga> */  -8.1142552237939071e+02);
  const double maxDispExpected =  (std::is_same_v<TestType, Grids::Yasp>) ? 0.786567027108460048 : ((std::is_same_v<TestType, Grids::Alu>)  ? 0.78426066482258983 :  /* std::is_same_v<TestType, Grids::Iga> */  0.615624125459537153);

  CHECK (energyExpected == Catch::Approx (nonLinOp.value()));
  CHECK (maxDispExpected == Catch::Approx (maxDisp));

  nonLinOp.template update<1>();
  CHECK (controlInfo.sucess);
  CHECK (gradTol >= nonLinOp.derivative().norm());
}