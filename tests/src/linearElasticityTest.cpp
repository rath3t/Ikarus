// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
using Dune::TestSuite;

#include "common.hh"
#include "easTest.hh"
#include "testHelpers.hh"

#include <numeric>
#include <string>
#include <vector>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/uggrid.hh>

#include <spdlog/spdlog.h>

#include <Eigen/Core>

#include "ikarus/utils/drawing/griddrawer.hh"
#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/mechanics/enhancedAssumedStrains.hh>
#include <ikarus/finiteElements/mechanics/linearElastic.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/utils/algorithms.hh>

/** This tests tests your element on some gridElement with some basis
 *
 * @tparam FEElementTemplate The element as template template parameter. The template needs to be the globalBasis
 * @tparam gridDim The dimension of the grid element the finite element should be tested
 * @tparam PreBasis The preBasis you want to test the element with
 * @tparam F A variadic number of the test functor you want to be checked, the need to accept a nonlinearoperator and
 * the finite element
 */
template <template <typename> typename FEElementTemplate, int gridDim, typename PreBasis, typename... F>
auto testFEElement(const PreBasis& preBasis, const std::string& elementName, F&&... f) {
  TestSuite t(std::string("testFEElement ") + elementName + " on grid element with dimension" + std::to_string(gridDim)
              + ".");

  auto fTuple = std::forward_as_tuple(f...);

  using Grid = Dune::UGGrid<gridDim>;

  std::vector<Dune::FieldVector<double, gridDim>> corners;

  const int numberOfVertices = Dune::power(2, gridDim);
  ValidCornerFactory<gridDim>::construct(corners, Dune::GeometryTypes::cube(gridDim));

  std::vector<unsigned int> vertexArrangment;
  vertexArrangment.resize(numberOfVertices);
  std::iota(vertexArrangment.begin(), vertexArrangment.end(), 0);

  Dune::GridFactory<Grid> gridFactory;
  for (auto& corner : corners) {
    gridFactory.insertVertex(corner);
  }
  gridFactory.insertElement(Dune::GeometryTypes::cube(gridDim), vertexArrangment);

  std::unique_ptr<Grid> grid = gridFactory.createGrid();

  auto gridView = grid->leafGridView();
  using namespace Ikarus;

  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, preBasis);

  auto localView = basis.localView();

  auto volumeLoad = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    VectorType fext;
    fext.setZero();
    fext[1] = 2 * lamb;
    fext[0] = lamb;
    return fext;
  };

  auto neumannBoundaryLoad = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    VectorType fext;
    fext.setZero();
    fext[1] = lamb / 40;
    fext[0] = 0;
    return fext;
  };

  // We artificially apply a Neumann load on the complete boundary
  Dune::BitSetVector<1> neumannVertices(gridView.size(2), true);

  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);

  const double youngsModulus = 1000;
  const double poissonsRatio = 0.3;
  auto element               = gridView.template begin<0>();

  using FEElementType = FEElementTemplate<decltype(basis)>;
  std::vector<FEElementType> fes;
  fes.emplace_back(basis, *element, youngsModulus, poissonsRatio, &volumeLoad, &neumannBoundary, &neumannBoundaryLoad);
  std::vector<bool> dirichletFlags(basis.size(), false);
  auto& fe             = fes[0];
  auto sparseAssembler = SparseFlatAssembler(basis, fes, dirichletFlags);

  typename FEElementType::FERequirementType::SolutionVectorType d;
  d.setRandom(basis.size());

  double lambda = 7.3;

  Ikarus::FErequirements requirements = FErequirementsBuilder()
                                            .insertParameter(Ikarus::FEParameter::loadfactor, lambda)
                                            .addAffordance(Ikarus::elastoStatics)
                                            .build();
  Eigen::VectorXd forces;
  Eigen::MatrixXd stiffnessmatrix;

  auto fvLambda = [&](auto&& d_) -> auto {
    forces.setZero(basis.localView().maxSize());
    requirements.setSolution(Ikarus::FESolutions::displacement, d_);
    return sparseAssembler.getScalar(requirements);
  };

  auto dfvLambda = [&](auto&& d_) -> auto& {
    forces.setZero(basis.localView().maxSize());
    requirements.setSolution(Ikarus::FESolutions::displacement, d_);
    return sparseAssembler.getVector(requirements);
  };
  auto ddfvLambda = [&](auto&& d_) -> auto& {
    requirements.setSolution(Ikarus::FESolutions::displacement, d_);
    return sparseAssembler.getMatrix(requirements);
  };
  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(fvLambda, dfvLambda, ddfvLambda), parameter(d));

  // execute all passed functions
  nonLinOp.updateAll();
  Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<sizeof...(F)>()),
                        [&](auto i) { t.subTest(std::get<i.value>(fTuple)(nonLinOp, fe)); });

  // check if element has a test functor, if yes we execute it
  if constexpr (requires { ElementTest<FEElementType>::test(); }) {
    auto testFunctor = ElementTest<FEElementType>::test();
    t.subTest(testFunctor(nonLinOp, fe));
  }

  return t;
}

template <typename Basis>
using EASElement = Ikarus::EnhancedAssumedStrains<Ikarus::LinearElastic<Basis>>;

template <typename Basis>
using LinearElasticElement = Ikarus::LinearElastic<Basis>;

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t("testFEElement");

  auto checkGradientFunctor
      = [](auto& nonLinOp, [[maybe_unused]] auto& fe) { return checkGradientOfElement(nonLinOp); };
  auto checkHessianFunctor  = [](auto& nonLinOp, [[maybe_unused]] auto& fe) { return checkHessianOfElement(nonLinOp); };
  auto checkJacobianFunctor = [](auto& nonLinOp, [[maybe_unused]] auto& fe) {
    auto subOperator = nonLinOp.template subOperator<1, 2>();
    return checkJacobianOfElement(subOperator);
  };
  using namespace Dune::Functions::BasisFactory;
  auto firstOrderLagrangePrePower2Basis  = power<2>(lagrange<1>(), FlatInterleaved());
  auto secondOrderLagrangePrePower2Basis = power<2>(lagrange<2>(), FlatInterleaved());
  auto firstOrderLagrangePrePower3Basis  = power<3>(lagrange<1>(), FlatInterleaved());
  auto secondOrderLagrangePrePower3Basis = power<3>(lagrange<2>(), FlatInterleaved());

  t.subTest(testFEElement<LinearElasticElement, 2>(firstOrderLagrangePrePower2Basis, "LinearElastic",
                                                   checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor));
  t.subTest(testFEElement<LinearElasticElement, 2>(secondOrderLagrangePrePower2Basis, "LinearElastic",
                                                   checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor));
  t.subTest(testFEElement<LinearElasticElement, 3>(firstOrderLagrangePrePower3Basis, "LinearElastic",
                                                   checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor));
  t.subTest(testFEElement<LinearElasticElement, 3>(secondOrderLagrangePrePower3Basis, "LinearElastic",
                                                   checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor));

  t.subTest(testFEElement<EASElement, 2>(firstOrderLagrangePrePower2Basis, "EAS", checkJacobianFunctor));
  t.subTest(testFEElement<EASElement, 3>(firstOrderLagrangePrePower3Basis, "EAS", checkJacobianFunctor));

  return t.exit();
}
