// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"
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

  Eigen::Matrix3d cc{
      { 0.600872, -0.179083, 0},
      {-0.179083,  0.859121, 0},
      {        0,         0, 1}
  };
  auto c = [&]() {
    if constexpr (strainTag == Ikarus::StrainTags::linear)
      return transformStrain<Ikarus::StrainTags::rightCauchyGreenTensor, Ikarus::StrainTags::greenLagrangian>(cc)
          .eval();
    else
      return transformStrain<Ikarus::StrainTags::rightCauchyGreenTensor, strainTag>(cc).eval();
  }();

  auto tol = isPlaneStress<IkarusMAT> ? 1e-8 : 1e-14;

  auto energy_muesli = muesliMat.template storedEnergy<strainTag>(c);
  auto stress_muesli = muesliMat.template stresses<strainTag>(c);
  auto moduli_muesli = muesliMat.template tangentModuli<strainTag>(c);

  auto energy_ikarus = ikarusMat.template storedEnergy<strainTag>(c);
  auto stress_ikarus = ikarusMat.template stresses<strainTag>(c);
  auto moduli_ikarus = ikarusMat.template tangentModuli<strainTag>(c);

  checkScalars(t, energy_muesli, energy_ikarus, "<energy<", tol);
  t.check(isApproxSame(stress_muesli, stress_ikarus, tol))
      << "Incorrect stresses." << " stress_muesli is\t" << stress_muesli.transpose() << "\n stress_ikarus is\t"
      << stress_ikarus.transpose();

  t.check(isApproxSame(moduli_muesli, moduli_ikarus, tol)) << "Incorrect tangentModuli." << " moduli_muesli is\n"
                                                           << moduli_muesli << "\n moduli_ikarus is\n"
                                                           << moduli_ikarus;

  return t;
}

template <typename MAT>
auto checkConstructors(Muesli::MaterialProperties matPar) {
  TestSuite t;
  auto mm = MAT(matPar);

  auto mm1{mm};
  auto mm2 = mm1;

  auto mm3 = std::forward<MAT>(mm);

  t.check(mm1.material().check()) << testLocation();
  t.check(mm2.material().check()) << testLocation();
  t.check(mm3.material().check()) << testLocation();

  t.check(mm1.assertMP()) << testLocation();
  t.check(mm2.assertMP()) << testLocation();
  t.check(mm3.assertMP()) << testLocation();

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  double Emod  = 1000;
  double nu    = 0.25; // Blatz Ko assumes nu = 0.25
  auto matPar_ = YoungsModulusAndPoissonsRatio{.emodul = Emod, .nu = nu};
  auto matPar  = toLamesFirstParameterAndShearModulus(matPar_);
  auto matProp = Muesli::propertiesFromIkarusMaterialParameters(matPar);

  auto lin  = LinearElasticity(matPar);
  auto linm = Muesli::SmallStrain(matPar);

  t.subTest(compareIkarusAndMuesli<StrainTags::linear>(linm, lin));

  auto nhm = Muesli::makeNeoHooke(matPar, false);
  auto nh  = NeoHooke(matPar);

  t.subTest(compareIkarusAndMuesli<StrainTags::rightCauchyGreenTensor>(nhm, nh));
  t.subTest(compareIkarusAndMuesli<StrainTags::deformationGradient>(nhm, nh));
  t.subTest(compareIkarusAndMuesli<StrainTags::greenLagrangian>(nhm, nh));

  auto svk  = StVenantKirchhoff(matPar);
  auto svkm = Muesli::makeSVK(matPar);

  t.subTest(compareIkarusAndMuesli<StrainTags::rightCauchyGreenTensor>(svkm, svk));
  t.subTest(compareIkarusAndMuesli<StrainTags::deformationGradient>(svkm, svk));
  t.subTest(compareIkarusAndMuesli<StrainTags::greenLagrangian>(svkm, svk));

  auto nhplaneStrain   = planeStrain(nh);
  auto nhplanseStrainm = planeStrain(nhm);
  t.subTest(compareIkarusAndMuesli<StrainTags::rightCauchyGreenTensor>(nhplanseStrainm, nhplaneStrain));

  auto nhplaneStress  = planeStress(nh);
  auto nhplaneStressm = planeStress(nhm);
  t.subTest(compareIkarusAndMuesli<StrainTags::rightCauchyGreenTensor>(nhplaneStressm, nhplaneStress));

  t.subTest(checkConstructors<Muesli::FiniteStrain<Muesli::NeoHooke>>(matProp));
  t.subTest(checkConstructors<Muesli::FiniteStrain<Muesli::StVenantKirchhoff>>(matProp));
  t.subTest(checkConstructors<Muesli::SmallStrain<Muesli::LinearElasticity>>(matProp));

  Muesli::makeNeoHooke(YoungsModulusAndBulkModulus{1000, 500});
  Muesli::makeNeoHooke(YoungsModulusAndPoissonsRatio{1000, 0.2});
  Muesli::makeSVK(YoungsModulusAndLamesFirstParameter{1000, 500});

  std::cout << svkm.name() << std::endl;
  std::cout << nhm.name() << std::endl;
  return t.exit();
}
