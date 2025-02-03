// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"
#include "testhyperelasticity.hh"
#include "tests/src/resultcollection.hh"

#include <muesli/muesli.h>

#include <dune/common/test/testsuite.hh>

#include <Eigen/Eigenvalues>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/materials/muesli/mueslimaterials.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/init.hh>

using namespace Ikarus;
using namespace Ikarus::Materials;
using Dune::TestSuite;

template <StrainTags strainTag, typename MuesliMAT, typename IkarusMAT>
auto compareIkarusAndMuesli(const MuesliMAT& muesliMat, const IkarusMAT& ikarusMat) {
  TestSuite t(MuesliMAT::name() + " vs " + IkarusMAT::name() + " InputStrainMeasure: " + toString(strainTag));

  Eigen::Matrix3d C = testMatrix();
  auto strain       = [&]() {
    if constexpr (strainTag == Ikarus::StrainTags::linear) // For linear we use the same as GreenLagrange
      return transformStrain<Ikarus::StrainTags::rightCauchyGreenTensor, Ikarus::StrainTags::greenLagrangian>(C).eval();
    else
      return transformStrain<Ikarus::StrainTags::rightCauchyGreenTensor, strainTag>(C).eval();
  }();

  auto tol = Testing::isPlaneStress<IkarusMAT> ? 1e-8 : 1e-14;

  auto energy_muesli = muesliMat.template storedEnergy<strainTag>(strain);
  auto stress_muesli = muesliMat.template stresses<strainTag>(strain);
  auto moduli_muesli = muesliMat.template tangentModuli<strainTag>(strain);

  auto energy_ikarus = ikarusMat.template storedEnergy<strainTag>(strain);
  auto stress_ikarus = ikarusMat.template stresses<strainTag>(strain);
  auto moduli_ikarus = ikarusMat.template tangentModuli<strainTag>(strain);

  checkScalars(t, energy_muesli, energy_ikarus, "Incorrect energy", tol);
  checkApproxMatrices(t, stress_muesli, stress_ikarus, "Incorrect stresses", tol);
  checkApproxMatrices(t, moduli_muesli, moduli_ikarus, "Incorrect tangentModuli", tol);

  return t;
}

template <typename MAT>
auto checkConstructors(muesli::materialProperties matPar) {
  TestSuite t;
  auto mm = MAT(matPar);

  auto mm1{mm};
  auto mm2 = mm1;

  auto mm3 = std::forward<MAT>(mm);

  t.check(mm1.material().check()) << testLocation();
  t.check(mm2.material().check()) << testLocation();
  t.check(mm3.material().check()) << testLocation();

  t.check(mm1.materialPoint().get() != NULL) << testLocation();
  t.check(mm2.materialPoint().get() != NULL) << testLocation();
  t.check(mm3.materialPoint().get() != NULL) << testLocation();

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  auto matPar_ = testMatPar();
  auto matPar  = toLamesFirstParameterAndShearModulus(matPar_);
  auto matProp = propertiesFromIkarusMaterialParameters(matPar);

  auto lin  = LinearElasticity(matPar);
  auto linm = SmallStrain(matPar);

  t.subTest(compareIkarusAndMuesli<StrainTags::linear>(linm, lin));

  auto nhm = makeMuesliNeoHooke(matPar, false);
  auto nh  = NeoHooke(matPar);

  t.subTest(compareIkarusAndMuesli<StrainTags::rightCauchyGreenTensor>(nhm, nh));
  t.subTest(compareIkarusAndMuesli<StrainTags::deformationGradient>(nhm, nh));
  t.subTest(compareIkarusAndMuesli<StrainTags::greenLagrangian>(nhm, nh));

  auto svk  = StVenantKirchhoff(matPar);
  auto svkm = makeMuesliSVK(matPar);

  t.subTest(compareIkarusAndMuesli<StrainTags::rightCauchyGreenTensor>(svkm, svk));
  t.subTest(compareIkarusAndMuesli<StrainTags::deformationGradient>(svkm, svk));
  t.subTest(compareIkarusAndMuesli<StrainTags::greenLagrangian>(svkm, svk));

  auto nhplaneStrain   = planeStrain(nh);
  auto nhplanseStrainm = planeStrain(nhm);
  t.subTest(compareIkarusAndMuesli<StrainTags::rightCauchyGreenTensor>(nhplanseStrainm, nhplaneStrain));

  auto nhplaneStress  = planeStress(nh);
  auto nhplaneStressm = planeStress(nhm);
  t.subTest(compareIkarusAndMuesli<StrainTags::rightCauchyGreenTensor>(nhplaneStressm, nhplaneStress));

  t.subTest(checkConstructors<FiniteStrain<muesli::neohookeanMaterial>>(matProp));
  t.subTest(checkConstructors<FiniteStrain<muesli::svkMaterial>>(matProp));
  t.subTest(checkConstructors<SmallStrain<muesli::elasticIsotropicMaterial>>(matProp));

  makeMuesliNeoHooke(YoungsModulusAndBulkModulus{1000, 500});
  makeMuesliNeoHooke(YoungsModulusAndPoissonsRatio{1000, 0.2});
  makeMuesliSVK(YoungsModulusAndLamesFirstParameter{1000, 500});

  t.check(svkm.name() == "FiniteStrain: SvkMaterial");
  t.check(nhm.name() == "FiniteStrain: NeohookeanMaterial");
  t.check(linm.name() == "SmallStrain: ElasticIsotropicMaterial");

  return t.exit();
}
