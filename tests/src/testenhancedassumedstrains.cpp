// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "resultcollection.hh"
#include "testeas.hh"
#include "testfeelement.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t("EAS");

  using namespace Dune::Functions::BasisFactory;
  auto firstOrderLagrangePrePower2Basis = power<2>(lagrange<1>(), FlatInterleaved());
  auto firstOrderLagrangePrePower3Basis = power<3>(lagrange<1>(), FlatInterleaved());
  constexpr auto randomlyDistorted      = CornerDistortionFlag::randomlyDistorted;
  constexpr auto unDistorted            = CornerDistortionFlag::unDistorted;

  auto linearElasticFunc3D = [](const Ikarus::YoungsModulusAndPoissonsRatio& parameter) {
    Ikarus::LinearElasticity lin(Ikarus::toLamesFirstParameterAndShearModulus(parameter));
    return Ikarus::linearElastic(lin);
  };
  auto linearElasticFuncPlaneStress = [](const Ikarus::YoungsModulusAndPoissonsRatio& parameter) {
    Ikarus::LinearElasticity lin(Ikarus::toLamesFirstParameterAndShearModulus(parameter));
    auto linPS = Ikarus::planeStress(lin, 1e-12);
    return Ikarus::linearElastic(linPS);
  };
  auto linearElasticFuncPlaneStrain = [](const Ikarus::YoungsModulusAndPoissonsRatio& parameter) {
    Ikarus::LinearElasticity lin(Ikarus::toLamesFirstParameterAndShearModulus(parameter));
    auto linPS = Ikarus::planeStrain(lin);
    return Ikarus::linearElastic(linPS);
  };

  t.subTest(testFEElement(firstOrderLagrangePrePower2Basis, "EAS", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::cube(), linearElasticFuncPlaneStress,
                          Ikarus::skills(Ikarus::eas()), Ikarus::AffordanceCollections::elastoStatics,
                          checkJacobianFunctor, singleElementTestFunctor));

  t.subTest(testFEElement(firstOrderLagrangePrePower3Basis, "EAS", randomlyDistorted,
                          Dune::ReferenceElements<double, 3>::cube(), linearElasticFunc3D,
                          Ikarus::skills(Ikarus::eas()), Ikarus::AffordanceCollections::elastoStatics,
                          checkJacobianFunctor, singleElementTestFunctor));

  // Plane Stress
  t.subTest(testFEElement(
      firstOrderLagrangePrePower2Basis, "EAS", unDistorted, Dune::ReferenceElements<double, 2>::cube(),
      linearElasticFuncPlaneStress, Ikarus::skills(Ikarus::eas()), Ikarus::AffordanceCollections::elastoStatics,
      checkCalculateAtFunctorFactory<Ikarus::ResultTypes::linearStress>(linearStressResultsOfSquare),
      checkCalculateAtFunctorFactory<Ikarus::ResultTypes::linearStress, false>(linearStressResultsOfSquare),
      checkResultFunctionFunctorFactory<Ikarus::ResultTypes::linearStress>(linearStressResultsOfSquare)));

  // Plane Strain
  t.subTest(testFEElement(
      firstOrderLagrangePrePower2Basis, "EAS", unDistorted, Dune::ReferenceElements<double, 2>::cube(),
      linearElasticFuncPlaneStrain, Ikarus::skills(Ikarus::eas()), Ikarus::AffordanceCollections::elastoStatics,
      checkCalculateAtFunctorFactory<Ikarus::ResultTypes::linearStress>(linearStressResultsOfSquare),
      checkCalculateAtFunctorFactory<Ikarus::ResultTypes::linearStress, false>(linearStressResultsOfSquare),
      checkResultFunctionFunctorFactory<Ikarus::ResultTypes::linearStress>(linearStressResultsOfSquare)));

  return t.exit();
}
