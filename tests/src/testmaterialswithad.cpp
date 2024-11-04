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
#include <ikarus/utils/init.hh>

using namespace Ikarus;
using Dune::TestSuite;

template <typename MAT, StrainTags strainTag>
auto stressByAD(const MAT& mat, const auto& c) {
  auto mat_ad = mat.template rebind<autodiff::dual>();

  auto f = [&](const auto& x) { return mat_ad.template storedEnergy<strainTag>(x); };

  auto dx = Eigen::Vector<autodiff::dual, 6>{};

  dx = toVoigt(c);
  autodiff::dual e;

  auto g = Eigen::Vector<double, 6>{};
  gradient(f, autodiff::wrt(dx), autodiff::at(dx), e, g);

  auto stress_nh_ad = (MAT::derivativeFactor * g).eval();

  return stress_nh_ad;
}

auto checkBlatzKo() {
  TestSuite t;

  // Eigen::Matrix3d e;
  // e.setRandom();
  // transformStrainAccordingToStrain<StrainTags::rightCauchyGreenTensor>(e);

  constexpr auto CauchyGreen = StrainTags::rightCauchyGreenTensor;

  // auto c = Eigen::Matrix3d::Identity().eval();
  Eigen::Matrix3d c{
      { 0.600872, -0.179083, 0},
      {-0.179083,  0.859121, 0},
      {        0,         0, 1}
  };

  // instantiate material models
  double Emod = 1000;
  double nu   = 0.25; // Blatz Ko assumes nu = 0.25
  auto matPar = YoungsModulusAndPoissonsRatio{.emodul = Emod, .nu = nu};
  auto mu     = convertLameConstants(matPar).toShearModulus();
  auto K      = convertLameConstants(matPar).toBulkModulus();
  auto Lambda = convertLameConstants(matPar).toLamesFirstParameter();

  auto nh = NeoHooke(toLamesFirstParameterAndShearModulus(matPar));

  auto energy_nh     = nh.storedEnergy<CauchyGreen>(c);
  auto stress_nh     = nh.stresses<CauchyGreen>(c);
  auto matTangent_nh = nh.tangentModuli<CauchyGreen>(c);

  std::cout << "Energy (NH):\n" << energy_nh << std::endl;
  std::cout << "Stress (NH):\n" << stress_nh << std::endl;
  // std::cout << "MatTangent (NH):\n" << matTangent_nh << std::endl;

  auto stress_nh_ad = stressByAD<decltype(nh), CauchyGreen>(nh, c);

  std::cout << "Stress (NH) AD:\n" << stress_nh_ad << std::endl;
  // std::cout << "MatTangent (NH):\n" << matTangent_nh << std::endl;

  std::array<double, 1> mu_og    = {mu};
  std::array<double, 1> alpha_og = {2.0};
  auto ogden_1                   = makeCompressibleOgden<1>(mu_og, alpha_og, {Lambda}, VF3{});

  auto stress_og_ad = stressByAD<decltype(ogden_1), CauchyGreen>(ogden_1, c);

  std::cout << "Stress (OG 1) AD:\n" << stress_og_ad << std::endl;

  auto ogden_2 = makeIncompressibleOgden<1>(mu_og, alpha_og, {Lambda}, VF3{});
  auto stress_og2 = ogden_2.stresses<CauchyGreen>(c);
  std::cout << "Stress (OG 2):\n" << stress_og2 << std::endl;

  auto stress_og2_ad = stressByAD<decltype(ogden_2), CauchyGreen>(ogden_2, c);
  std::cout << "Stress (OG 2) AD :\n" << stress_og2_ad << std::endl;

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  LamesFirstParameterAndShearModulus matPar{.lambda = 1000, .mu = 500};

  t.subTest(checkBlatzKo());

  return t.exit();
}
