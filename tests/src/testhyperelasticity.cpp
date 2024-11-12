// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Eigenvalues>

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <autodiff/forward/utils/derivative.hpp>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/functionsanitychecks.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/nonlinearoperator.hh>

using namespace Ikarus;
using namespace Materials;
using Dune::TestSuite;

inline auto testMatrix() {
  return Eigen::Matrix3d{
      { 0.600872, -0.179083, 0},
      {-0.179083,  0.859121, 0},
      {        0,         0, 1}
  };
}

auto testMatPar() {
  double Emod = 1000;
  double nu   = 0.25; // Blatz Ko assumes nu = 0.25
  return YoungsModulusAndPoissonsRatio{.emodul = Emod, .nu = nu};
}

auto testVolumetricFunctions() {
  const auto detC = testMatrix().determinant();
  const auto J    = sqrt(detC);

  TestSuite t("Test volumetric functions by AD");

  auto checkFirstAndSecondDerivativeOfEnergy = [&]<typename VF>(const VF& vf) {
    auto vf_ad          = vf.template rebind<autodiff::dual2nd>();
    autodiff::dual2nd x = J;
    auto f              = [&](const auto& xloc) { return vf_ad.storedEnergyImpl(xloc); };

    auto derivs = derivatives(f, autodiff::wrt(x), autodiff::at(x));

    auto uprime      = vf.firstDerivativeImpl(J);
    auto uprimeprime = vf.secondDerivativeImpl(J);

    t.check(Dune::FloatCmp::eq(uprime, derivs[1], 1e-14))
        << testLocation() << "Uprime is " << uprime << " but should be " << derivs[1] << " according to AD for "
        << Dune::className<VF>();
    t.check(Dune::FloatCmp::eq(uprimeprime, derivs[2], 1e-14))
        << testLocation() << "Uprimeprime is " << uprimeprime << " but should be " << derivs[2]
        << " according to AD for " << Dune::className<VF>();
  };

  auto checkFirstDerivativeOfFirstDerivative = [&]<typename VF>(const VF& vf) {
    auto vf_ad       = vf.template rebind<autodiff::dual>();
    autodiff::dual x = J;
    auto f           = [&](const auto& xloc) { return vf_ad.firstDerivativeImpl(xloc); };

    auto uprimeprime_ad = derivative(f, autodiff::wrt(x), autodiff::at(x));

    auto uprimeprime = vf.secondDerivativeImpl(J);

    t.check(Dune::FloatCmp::eq(uprimeprime, uprimeprime_ad, 1e-14))
        << testLocation() << "Uprimeprime is " << uprimeprime << " but should be " << uprimeprime_ad
        << " according to AD for " << Dune::className<VF>() << " by taking first derivative of first derivative of U";
  };

  auto vf1 = VF1{};
  checkFirstAndSecondDerivativeOfEnergy(vf1);
  checkFirstDerivativeOfFirstDerivative(vf1);

  auto vf2 = VF2{};
  checkFirstAndSecondDerivativeOfEnergy(vf2);
  checkFirstDerivativeOfFirstDerivative(vf2);

  auto vf3 = VF3{};
  checkFirstAndSecondDerivativeOfEnergy(vf3);
  checkFirstDerivativeOfFirstDerivative(vf3);

  auto vf4 = VF4{0.5};
  checkFirstAndSecondDerivativeOfEnergy(vf4);
  checkFirstDerivativeOfFirstDerivative(vf4);

  auto vf5 = VF5{};
  checkFirstAndSecondDerivativeOfEnergy(vf5);
  checkFirstDerivativeOfFirstDerivative(vf5);

  auto vf6 = VF6{};
  checkFirstAndSecondDerivativeOfEnergy(vf6);
  checkFirstDerivativeOfFirstDerivative(vf6);

  return t;
}

auto recoverNeoHookeThroughOgden() {
  using StrainTags::rightCauchyGreenTensor;

  TestSuite t("Recover NeoHooke material through Ogden test");

  auto matPar = testMatPar();
  auto mu     = convertLameConstants(matPar).toShearModulus();
  auto Lambda = convertLameConstants(matPar).toLamesFirstParameter();
  auto c      = testMatrix();

  std::array<double, 1> mu_og    = {mu};
  std::array<double, 1> alpha_og = {2.0};

  auto ogden = makeOgden<1, PrincipalStretchTag::total>(mu_og, alpha_og, {Lambda}, VF3{});
  auto nh    = NeoHooke(toLamesFirstParameterAndShearModulus(matPar));

  auto energy_og = ogden.storedEnergy<rightCauchyGreenTensor>(c);
  auto stress_og = ogden.stresses<rightCauchyGreenTensor>(c);
  auto moduli_og = ogden.tangentModuli<rightCauchyGreenTensor>(c);

  auto energy_nh = nh.storedEnergy<rightCauchyGreenTensor>(c);
  auto stress_nh = nh.stresses<rightCauchyGreenTensor>(c);
  auto moduli_nh = nh.tangentModuli<rightCauchyGreenTensor>(c);

  t.check(Dune::FloatCmp::eq(energy_og, energy_nh, 1e-14)) << testLocation() << "Energy not the same";
  t.check(isApproxSame(stress_og, stress_nh, 1e-14)) << testLocation() << "Stress not the same";
  t.check(isApproxSame(moduli_og, moduli_nh, 1e-14)) << testLocation() << "Tangent not the same";

  return t;
}

// auto playground() {
//   TestSuite t;

// // Eigen::Matrix3d e;
// // e.setRandom();
// // transformStrainAccordingToStrain<StrainTags::rightCauchyGreenTensor>(e);

// constexpr auto CauchyGreen = StrainTags::rightCauchyGreenTensor;

// // auto c = Eigen::Matrix3d::Identity().eval();
// auto c = testMatrix();

// // instantiate material models
// double Emod = 1000;
// double nu   = 0.25; // Blatz Ko assumes nu = 0.25
// auto matPar = YoungsModulusAndPoissonsRatio{.emodul = Emod, .nu = nu};
// auto mu     = convertLameConstants(matPar).toShearModulus();
// auto K      = convertLameConstants(matPar).toBulkModulus();
// auto Lambda = convertLameConstants(matPar).toLamesFirstParameter();

// // auto nh = NeoHooke(toLamesFirstParameterAndShearModulus(matPar));

// // auto energy_nh     = nh.storedEnergy<CauchyGreen>(c);
// // auto stress_nh     = nh.stresses<CauchyGreen>(c);
// // auto matTangent_nh = nh.tangentModuli<CauchyGreen>(c);

// // std::cout << "Energy (NH):\n" << energy_nh << std::endl;
// // std::cout << "Stress (NH):\n" << stress_nh << std::endl;
// // std::cout << "MatTangent (NH):\n" << matTangent_nh << std::endl;

// auto bk = makeBlatzKo(ShearModulus{mu});

// // auto energy_bk     = bk.storedEnergy<CauchyGreen>(c);
// // auto stress_bk     = bk.stresses<CauchyGreen>(c);
// // auto matTangent_bk = bk.tangentModuli<CauchyGreen>(c);

// // std::cout << "Energy (BK):\n" << energy_bk << std::endl;
// // std::cout << "Stress (BK):\n" << stress_bk << std::endl;
// // std::cout << "MatTangent (BK):\n" << matTangent_bk << std::endl;

// std::array<double, 1> mu_og    = {mu};
// std::array<double, 1> alpha_og = {2.0};

// auto ogden_1    = makeModifiedOgden<1>(mu_og, alpha_og, {Lambda}, VF3{});
// auto energy_og1 = ogden_1.storedEnergy<CauchyGreen>(c);
// auto stress_og1 = ogden_1.stresses<CauchyGreen>(c);
// auto moduli_og1 = ogden_1.tangentModuli<CauchyGreen>(c);

// std::cout << "Energy (OG 1):\n" << energy_og1 << std::endl;
// std::cout << "Stress (OG 1):\n" << stress_og1 << std::endl;
// std::cout << "MatTangent (OG 1):\n" << moduli_og1 << std::endl;

// // test with voigt notation
// auto cvoigt      = toVoigt(c, true);
// auto stress_og1v = ogden_1.stresses<CauchyGreen>(cvoigt);

// assert(isApproxSame(stress_og1, stress_og1v, 1e-10));

// auto ogden_2 = makeRegularizedOgden<1>(mu_og, alpha_og, {K}, VF3{});

// auto energy_og2 = ogden_2.storedEnergy<CauchyGreen>(c);
// auto stress_og2 = ogden_2.stresses<CauchyGreen>(c);
// auto moduli_og2 = ogden_2.tangentModuli<CauchyGreen>(c);

// std::cout << "Energy (OG 2):\n" << energy_og2 << std::endl;
// std::cout << "Stress (OG 2):\n" << stress_og2 << std::endl;
// std::cout << "MatTangent (OG 2):\n" << moduli_og2 << std::endl;

// auto ogden_3 = makeRegularizedOgden<1>(mu_og, alpha_og);

// std::cout << bk.name() << std::endl;
// std::cout << ogden_1.name() << std::endl;
// std::cout << ogden_2.name() << std::endl;
// std::cout << ogden_3.name() << std::endl;

// using Ikarus::RegularizedTag;
// auto ogden_4 = makeOgden<1, RegularizedTag::regularized>(mu_og, alpha_og);
// static_assert(std::same_as<decltype(ogden_4), decltype(ogden_3)>, "Should be same");

// auto ogden_5 = makeOgden<1, RegularizedTag::modified>(mu_og, alpha_og, {Lambda}, VF3{});
// static_assert(std::same_as<decltype(ogden_5), decltype(ogden_1)>, "Should be same");

// auto psOgden  = planeStrain(ogden_5);
// auto pstOgden = planeStress(ogden_5);

// auto stress_og5ps  = psOgden.stresses<CauchyGreen>(c);
// auto stress_og5pst = pstOgden.stresses<CauchyGreen>(c);

// std::cout << "Stress (OG PS):\n" << stress_og5ps << std::endl;
// std::cout << "Stress (OG PST):\n" << stress_og5pst << std::endl;

// return t;
// }

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  // t.subTest(playground());

  t.subTest(testVolumetricFunctions());
  t.subTest(recoverNeoHookeThroughOgden());

  return t.exit();
}
