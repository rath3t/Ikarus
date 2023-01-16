// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

using Dune::TestSuite;

#include "testFEElement.hh"

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <ikarus/finiteElements/mechanics/linearElastic.hh>
#include <ikarus/utils/duneUtilities.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>

template <typename Basis>
using LinearElasticElement = Ikarus::LinearElastic<Basis>;

auto checkCauchyStress2D() {
  TestSuite t("Cauchy Stress check for the 2D case");
  using namespace Ikarus;
  using namespace Dune::Indices;
  using namespace Dune::Functions::BasisFactory;
  constexpr int gridDim = 2;

  using Grid = Dune::YaspGrid<gridDim>;

  Dune::FieldVector<double, 2> bbox       = {1, 1};
  std::array<int, 2> elementsPerDirection = {1, 1};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);
  auto gridView                           = grid->leafGridView();
  auto basis = Ikarus::makeConstSharedBasis(gridView, power<gridDim>(lagrange<1>(), FlatInterleaved()));
  std::vector<Ikarus::LinearElastic<typename decltype(basis)::element_type>> fes;

  for (auto& element : elements(gridView)) {
    auto localView = basis->localView();
    fes.emplace_back(*basis, element, 100, 0.0);
  }

  Eigen::VectorXd displacement;
  displacement.resize(basis->size());
  displacement << 0, 0, 1, 1, 1, 1, 1, 1;
  std::vector<Dune::FieldVector<double, 3>> expectedStress = {{100, 100, 100}, {100, 0, 50}, {0, 100, 50}, {0, 0, 0}};

  auto feRequirements = FErequirements().addAffordance(Ikarus::AffordanceCollections::elastoStatics);
  feRequirements.insertGlobalSolution(Ikarus::FESolutions::displacement, displacement);

  auto resultRequirements = Ikarus::ResultRequirements<Eigen::VectorXd>()
                                .insertGlobalSolution(Ikarus::FESolutions::displacement, displacement)
                                .addResultRequest(ResultType::cauchyStress);
  ResultTypeMap<double> result;

  auto scalarBasis        = makeConstSharedBasis(gridView, lagrangeDG<1>());
  auto localScalarView    = scalarBasis->localView();
  bool checkStressAtNodes = false;
  std::vector<Dune::FieldVector<double, 3>> stressVector(scalarBasis->size());
  auto ele = elements(gridView).begin();

  for (auto& fe : fes) {
    localScalarView.bind(*ele);
    const auto& fe2              = localScalarView.tree().finiteElement();
    const auto& referenceElement = Dune::ReferenceElements<double, gridDim>::general(ele->type());
    for (auto c = 0UL; c < fe2.size(); ++c) {
      const auto fineKey                        = fe2.localCoefficients().localKey(c);
      const auto nodalPositionInChildCoordinate = referenceElement.position(fineKey.subEntity(), fineKey.codim());
      auto coord                                = toEigen(nodalPositionInChildCoordinate);
      fe.calculateAt(resultRequirements, coord, result);
      Eigen::Vector3d computedResult = result.getResult(ResultType::cauchyStress).eval();
      stressVector[localScalarView.index(localScalarView.tree().localIndex(c))[0]] = toDune(computedResult);
      for (auto voigtSize = 0UL; voigtSize < 3; ++voigtSize)
        if (Dune::FloatCmp::eq(stressVector[localScalarView.index(localScalarView.tree().localIndex(c))[0]][voigtSize],
                               expectedStress[c][voigtSize]))
          checkStressAtNodes = true;
        else {
          checkStressAtNodes = false;
          break;
        }
      if (not checkStressAtNodes) break;
    }
    ++ele;
  }
  t.check(checkStressAtNodes, "Stresses at nodes are correctly computed");
  return t;
}

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t("LinearElasticity");

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
  t.subTest(checkCauchyStress2D());

  return t.exit();
}
