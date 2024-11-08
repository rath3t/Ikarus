// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "resultcollection.hh"
#include "testfeelement.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <ikarus/finiteelements/feresulttypes.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mechanics/materials/muesli/mueslimaterials.hh>
#include <ikarus/io/resultevaluators.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

int main(int argc, char** argv) {
  using namespace Ikarus;
  using namespace ResultTypes;
  using namespace ResultEvaluators;

  init(argc, argv);
  TestSuite t("LinearElasticity");

  using namespace Dune::Functions::BasisFactory;
  using namespace Ikarus::Materials;

  auto firstOrderLagrangePrePower2Basis         = power<2>(lagrange<1>(), FlatInterleaved());
  auto secondOrderLagrangePrePower2Basis        = power<2>(lagrange<2>(), FlatInterleaved());
  auto firstOrderLagrangePrePower3Basis         = power<3>(lagrange<1>(), FlatInterleaved());
  auto secondOrderLagrangePrePower3Basis        = power<3>(lagrange<2>(), FlatInterleaved());
  auto secondOrderLagrangePrePower3BasisBlocked = power<3>(lagrange<2>());
  constexpr auto randomlyDistorted              = CornerDistortionFlag::randomlyDistorted;
  constexpr auto unDistorted                    = CornerDistortionFlag::unDistorted;

  // Test cube 2D
  auto linearElasticFunc3D = [](const Ikarus::YoungsModulusAndPoissonsRatio& parameter) {
    LinearElasticity lin(Ikarus::toLamesFirstParameterAndShearModulus(parameter));
    return Ikarus::linearElastic(lin);
  };
  auto linearElasticFunc3D_Muesli = [](const Ikarus::YoungsModulusAndPoissonsRatio& parameter) {
    auto lin = Ikarus::Materials::Muesli::makeLinearElasticity(parameter);
    return Ikarus::linearElastic(lin);
  };
  auto linearElasticFuncPlaneStress = [](const Ikarus::YoungsModulusAndPoissonsRatio& parameter) {
    LinearElasticity lin(Ikarus::toLamesFirstParameterAndShearModulus(parameter));
    auto linPS = planeStress(lin);
    return Ikarus::linearElastic(linPS);
  };
  auto linearElasticFuncPlaneStrain = [](const Ikarus::YoungsModulusAndPoissonsRatio& parameter) {
    LinearElasticity lin(Ikarus::toLamesFirstParameterAndShearModulus(parameter));
    auto linPS = planeStrain(lin);
    return Ikarus::linearElastic(linPS);
  };

  // Plane stress
  t.subTest(testFEElement(firstOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::cube(), linearElasticFuncPlaneStress, skills(),
                          AffordanceCollections::elastoStatics, checkGradientFunctor, checkHessianFunctor,
                          checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement(secondOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::cube(), linearElasticFuncPlaneStress, skills(),
                          AffordanceCollections::elastoStatics, checkGradientFunctor, checkHessianFunctor,
                          checkJacobianFunctor, checkFEByAutoDiffFunctor));

  t.subTest(testFEElement(
      firstOrderLagrangePrePower2Basis, "LinearElastic", unDistorted, Dune::ReferenceElements<double, 2>::cube(),
      linearElasticFuncPlaneStress, skills(), AffordanceCollections::elastoStatics,
      checkCalculateAtFunctorFactory<linearStress>(linearStressResultsOfSquare),
      checkCalculateAtFunctorFactory<linearStress, false>(linearStressResultsOfSquare),
      checkResultFunctionFunctorFactory<linearStress>(linearStressResultsOfSquare),
      checkResultFunctionFunctorFactory<linearStress, PolarStress>(linearPolarStressResultsOfSquare,
                                                                   Dune::FieldVector<double, 2>{0.5, 0.5}),
      checkResultFunctionFunctorFactory<linearStress, VonMises>(linearVonMisesResultsOfSquare),
      checkResultFunctionFunctorFactory<linearStress, HydrostaticStress>(linearHydrostaticStressResultsOfSquare),
      checkResultFunctionFunctorFactory<linearStress, Triaxiality>(linearTriaxialityStressResultsOfSquare),
      checkResultFunctionFunctorFactory<linearStress, PrincipalStress<2>>(linearPrincipalStressResultsOfSquare)));

  // Plane strain
  t.subTest(testFEElement(firstOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::cube(), linearElasticFuncPlaneStrain, skills(),
                          AffordanceCollections::elastoStatics, checkGradientFunctor, checkHessianFunctor,
                          checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement(secondOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::cube(), linearElasticFuncPlaneStrain, skills(),
                          AffordanceCollections::elastoStatics, checkGradientFunctor, checkHessianFunctor,
                          checkJacobianFunctor, checkFEByAutoDiffFunctor));

  t.subTest(testFEElement(
      firstOrderLagrangePrePower2Basis, "LinearElastic", unDistorted, Dune::ReferenceElements<double, 2>::cube(),
      linearElasticFuncPlaneStrain, skills(), AffordanceCollections::elastoStatics,
      checkCalculateAtFunctorFactory<linearStress>(linearStressResultsOfSquare),
      checkCalculateAtFunctorFactory<linearStressFull>(linear3dPlaneStrainStressResultsOfSquare),
      checkCalculateAtFunctorFactory<linearStress, false>(linearStressResultsOfSquare),
      // checkCalculateAtFunctorFactory<linearStressFull, false>(linear3dPlaneStrainStressResultsOfSquare),
      checkResultFunctionFunctorFactory<linearStress>(linearStressResultsOfSquare),
      checkResultFunctionFunctorFactory<linearStressFull>(linear3dPlaneStrainStressResultsOfSquare),
      checkResultFunctionFunctorFactory<linearStressFull, VonMises>(linearVonMisesResultsOfSquare),
      checkResultFunctionFunctorFactory<linearStressFull, HydrostaticStress>(linearHydrostaticStressResultsOfSquare),
      checkResultFunctionFunctorFactory<linearStressFull, Triaxiality>(linearTriaxialityStressResultsOfSquare),
      checkResultFunctionFunctorFactory<linearStressFull, PrincipalStress<3>>(linearPrincipalStressResultsOfSquare)));

  // Test simplex 2D
  t.subTest(testFEElement(firstOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::simplex(), linearElasticFuncPlaneStress, skills(),
                          AffordanceCollections::elastoStatics, checkGradientFunctor, checkHessianFunctor,
                          checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement(secondOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::simplex(), linearElasticFuncPlaneStress, skills(),
                          AffordanceCollections::elastoStatics, checkGradientFunctor, checkHessianFunctor,
                          checkJacobianFunctor, checkFEByAutoDiffFunctor));

  t.subTest(testFEElement(firstOrderLagrangePrePower2Basis, "LinearElastic", unDistorted,
                          Dune::ReferenceElements<double, 2>::simplex(), linearElasticFuncPlaneStress, skills(),
                          AffordanceCollections::elastoStatics,
                          checkCalculateAtFunctorFactory<linearStress>(linearStressResultsOfTriangle),
                          checkResultFunctionFunctorFactory<linearStress, PolarStress>(
                              linearPolarStressResultsOfTriangle, Dune::FieldVector<double, 2>{1.0 / 3, 1.0 / 3})));

  t.subTest(testFEElement(firstOrderLagrangePrePower2Basis, "LinearElastic", unDistorted,
                          Dune::ReferenceElements<double, 2>::simplex(), linearElasticFuncPlaneStrain, skills(),
                          AffordanceCollections::elastoStatics,
                          checkCalculateAtFunctorFactory<linearStress>(linearStressResultsOfTriangle)));

  // Test cube 3D
  t.subTest(testFEElement(firstOrderLagrangePrePower3Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 3>::cube(), linearElasticFunc3D, skills(),
                          AffordanceCollections::elastoStatics, checkGradientFunctor, checkHessianFunctor,
                          checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement(secondOrderLagrangePrePower3Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 3>::cube(), linearElasticFunc3D, skills(),
                          AffordanceCollections::elastoStatics, checkGradientFunctor, checkHessianFunctor,
                          checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement(secondOrderLagrangePrePower3BasisBlocked, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 3>::cube(), linearElasticFunc3D, skills(),
                          AffordanceCollections::elastoStatics, checkGradientFunctor, checkHessianFunctor,
                          checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement(firstOrderLagrangePrePower3Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 3>::cube(), linearElasticFunc3D_Muesli, Ikarus::skills(),
                          Ikarus::AffordanceCollections::elastoStatics, checkGradientFunctor, checkHessianFunctor,
                          checkJacobianFunctor));

  t.subTest(testFEElement(
      firstOrderLagrangePrePower3Basis, "LinearElastic", unDistorted, Dune::ReferenceElements<double, 3>::cube(),
      linearElasticFunc3D, skills(), AffordanceCollections::elastoStatics,
      checkCalculateAtFunctorFactory<linearStress>(linearStressResultsOfCube),
      checkResultFunctionFunctorFactory<linearStress>(linearStressResultsOfCube),
      checkResultFunctionFunctorFactory<linearStress, VonMises>(linearVonMisesResultsOfCube),
      checkResultFunctionFunctorFactory<linearStress, HydrostaticStress>(linearHydrostaticStressResultsOfCube),
      checkResultFunctionFunctorFactory<linearStress, Triaxiality>(linearTriaxialityResultsOfCube),
      checkResultFunctionFunctorFactory<linearStress, PrincipalStress<3>>(linearPrincipalStressResultsOfCube)));

  // Test simplex 3D
  t.subTest(testFEElement(firstOrderLagrangePrePower3Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 3>::simplex(), linearElasticFunc3D, skills(),
                          AffordanceCollections::elastoStatics, checkGradientFunctor, checkHessianFunctor,
                          checkJacobianFunctor, checkFEByAutoDiffFunctor));

  t.subTest(testFEElement(firstOrderLagrangePrePower3Basis, "LinearElastic", unDistorted,
                          Dune::ReferenceElements<double, 3>::simplex(), linearElasticFunc3D, skills(),
                          AffordanceCollections::elastoStatics,
                          checkCalculateAtFunctorFactory<linearStress>(linearStressResultsOfTetrahedron),
                          checkResultFunctionFunctorFactory<linearStress>(linearStressResultsOfTetrahedron)));
  return t.exit();
}
