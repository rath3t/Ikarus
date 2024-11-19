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
  const auto detC      = testMatrix().determinant();
  const auto J         = sqrt(detC);
  constexpr double tol = 1e-14;

  TestSuite t("Test volumetric functions by AD");

  auto checkFirstAndSecondDerivativeOfEnergy = [&]<typename VF>(const VF& vf) {
    auto vf_ad          = vf.template rebind<autodiff::dual2nd>();
    autodiff::dual2nd x = J;
    auto f              = [&](const auto& xloc) { return vf_ad.storedEnergyImpl(xloc); };
    auto derivs         = derivatives(f, autodiff::wrt(x), autodiff::at(x));
    auto uprime         = vf.firstDerivativeImpl(J);
    auto uprimeprime    = vf.secondDerivativeImpl(J);

    checkScalars(t, uprime, derivs[1], testLocation() + Dune::className<VF>(), tol);
    checkScalars(t, uprimeprime, derivs[2], testLocation() + Dune::className<VF>(), tol);
  };

  auto checkFirstDerivativeOfFirstDerivative = [&]<typename VF>(const VF& vf) {
    auto vf_ad          = vf.template rebind<autodiff::dual>();
    autodiff::dual x    = J;
    auto f              = [&](const auto& xloc) { return vf_ad.firstDerivativeImpl(xloc); };
    auto uprimeprime_ad = derivative(f, autodiff::wrt(x), autodiff::at(x));
    auto uprimeprime    = vf.secondDerivativeImpl(J);

    checkScalars(t, uprimeprime, uprimeprime_ad, testLocation() + Dune::className<VF>(), tol);
  };

  auto checkVolumetricFunctionsByAD = [&]<typename VF>(const VF& vf) {
    checkFirstAndSecondDerivativeOfEnergy(vf);
    checkFirstDerivativeOfFirstDerivative(vf);
  };

  checkVolumetricFunctionsByAD(VF1{});
  checkVolumetricFunctionsByAD(VF2{});
  checkVolumetricFunctionsByAD(VF3{});
  checkVolumetricFunctionsByAD(VF4{0.5});
  checkVolumetricFunctionsByAD(VF5{});
  checkVolumetricFunctionsByAD(VF6{});
  checkVolumetricFunctionsByAD(VF7{0.5});
  checkVolumetricFunctionsByAD(VF8{});
  checkVolumetricFunctionsByAD(VF9{});
  checkVolumetricFunctionsByAD(VF10{0.4});
  checkVolumetricFunctionsByAD(VF11{});

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

  constexpr double tol = 1e-14;

  checkScalars(t, energy_og, energy_nh, testLocation() + "Incorrect Energy.", tol);
  t.check(isApproxSame(stress_og, stress_nh, tol))
      << testLocation() << "Incorrect stresses." << " stress_og is\t" << stress_og.transpose() << "\n stress_nh is\t"
      << stress_nh.transpose();
  t.check(isApproxSame(moduli_og, moduli_nh, tol)) << testLocation() << "Incorrect tangentModuli." << " moduli_og is\n"
                                                   << moduli_og << "\n moduli_nh is\n"
                                                   << moduli_nh;

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testVolumetricFunctions());
  t.subTest(recoverNeoHookeThroughOgden());

  return t.exit();
}
