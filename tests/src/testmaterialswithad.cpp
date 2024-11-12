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
using namespace Ikarus::Materials;
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

  auto stress_ad = (MAT::derivativeFactor * g).eval();

  return stress_ad;
}

template <typename MAT, StrainTags strainTag>
auto mattangentByAD(const MAT& mat, const auto& c) {
  auto mat_ad = mat.template rebind<autodiff::dual>();

  auto f = [&](const auto& x) { return mat_ad.template stresses<strainTag>(x); };

  auto dx = Eigen::Vector<autodiff::dual, 6>{};

  dx = toVoigt(c);
  Eigen::VectorXdual g(6);

  auto h = Eigen::Matrix<double, 6, 6>{};
  jacobian(f, autodiff::wrt(dx), autodiff::at(dx), g, h);

  auto matTangent_ad = (MAT::derivativeFactor * h).eval();

  return matTangent_ad;
}

template <typename MAT, StrainTags strainTag>
auto mattangentByADWithEnergy(const MAT& mat, const auto& c) {
  auto mat_ad = mat.template rebind<autodiff::dual2nd>();

  auto f = [&](const auto& x) { return mat_ad.template storedEnergy<strainTag>(x); };

  Eigen::Matrix<autodiff::dual2nd, 6, 1> dx = toVoigt(c);

  autodiff::dual2nd e;
  Eigen::Matrix<double, 6, 1> g;
  Eigen::Matrix<double, 6, 6> h;

  h = autodiff::hessian(f, autodiff::wrt(dx), autodiff::at(dx), e, g);

  auto matTangent_ad = (MAT::derivativeFactor * MAT::derivativeFactor * h).eval();

  return matTangent_ad;
}

template <StrainTags straintag, typename MAT>
auto testMaterial(const MAT& mat, const auto& c, double prec = 1e-8) {
  TestSuite t("Test Material with AD: " + MAT::name() /*, Dune::TestSuite::AlwaysThrow */);

  auto stress     = mat.template stresses<straintag>(c);
  auto matTangent = mat.template tangentModuli<straintag>(c);

  auto stress_ad       = stressByAD<MAT, straintag>(mat, c);
  auto matTangent_ad   = mattangentByAD<MAT, straintag>(mat, c);
  auto matTangent_ad_e = mattangentByADWithEnergy<MAT, straintag>(mat, c);

  t.check(isApproxSame(stress, stress_ad, prec)) << "stress wrong";
  t.check(isApproxSame(matTangent, matTangent_ad, prec)) << "matTangent wrong";
  t.check(isApproxSame(matTangent, matTangent_ad_e, prec)) << "matTangent from stress wrong";
  t.check(isApproxSame(matTangent_ad_e, matTangent_ad, prec))
      << "matTangent from stress not equal to mattangent from energy";
  ;

  return t;
}

// auto testMaterialByAD() {
//   TestSuite t;

// Eigen::Matrix3d c{
//     { 0.600872, -0.179083, 0},
//     {-0.179083,  0.859121, 0},
//     {        0,         0, 1}
// };

// constexpr auto CauchyGreen = StrainTags::rightCauchyGreenTensor;

// // instantiate material models
// double Emod = 1000;
// double nu   = 0.25; // Blatz Ko assumes nu = 0.25
// auto matPar = YoungsModulusAndPoissonsRatio{.emodul = Emod, .nu = nu};
// auto mu     = convertLameConstants(matPar).toShearModulus();
// auto K      = convertLameConstants(matPar).toBulkModulus();
// auto Lambda = convertLameConstants(matPar).toLamesFirstParameter();

// auto nh = NeoHooke(toLamesFirstParameterAndShearModulus(matPar));

// // auto stress_nh_ad       = stressByAD<decltype(nh), CauchyGreen>(nh, c);
// // auto matTangent_nh_ad   = mattangentByAD<decltype(nh), CauchyGreen>(nh, c);
// // auto matTangent_nh_ad_e = mattangentByADWithEnergy<decltype(nh), CauchyGreen>(nh, c);

// // std::cout << "Stress (NH) AD:\n" << stress_nh_ad << std::endl;
// // std::cout << "MatTangent (NH) AD:\n" << matTangent_nh_ad << std::endl;
// // std::cout << "MatTangent (NH) AD Energy:\n" << matTangent_nh_ad_e << std::endl;

// std::array<double, 1> mu_og    = {mu};
// std::array<double, 1> alpha_og = {2.0};
// // auto ogden_1                   = makeModifiedOgden<1>(mu_og, alpha_og, {Lambda}, VF3{});

// // auto energy_og     = ogden_1.storedEnergy<CauchyGreen>(c);
// // auto stress_og     = ogden_1.stresses<CauchyGreen>(c);
// // auto matTangent_og = ogden_1.tangentModuli<CauchyGreen>(c);

// // auto stress_og1_ad        = stressByAD<decltype(ogden_1), CauchyGreen>(ogden_1, c);
// // auto matTangent_og1_ad   = mattangentByAD<decltype(ogden_1), CauchyGreen>(ogden_1, c);
// // auto matTangent_og1_ad_e = mattangentByADWithEnergy<decltype(ogden_1), CauchyGreen>(ogden_1, c);

// auto ogden_2 = makeOgden<1, PrincipalStretchTag::total>(mu_og, alpha_og, {Lambda}, VF3{});

// auto stress_og2    = ogden_2.stresses<CauchyGreen>(c);
// auto matTangent_og = ogden_2.tangentModuli<CauchyGreen>(c);

// std::cout << "Stress (OG 2):\n" << stress_og2 << std::endl;
// std::cout << "MatTangent (OG 2):\n" << matTangent_og << std::endl;

// auto stress_og2_ad       = stressByAD<decltype(ogden_2), CauchyGreen>(ogden_2, c);
// auto matTangent_og2_ad_e = mattangentByADWithEnergy<decltype(ogden_2), CauchyGreen>(ogden_2, c);
// auto matTangent_og2_ad   = mattangentByAD<decltype(ogden_2), CauchyGreen>(ogden_2, c);

// std::cout << "Stress (OG 2) AD:\n" << stress_og2_ad << std::endl;
// std::cout << "MatTangent (OG 2) AD via energy:\n" << matTangent_og2_ad_e << std::endl;
// // std::cout << "MatTangent (OG 2) AD via stresses:\n" << matTangent_og2_ad << std::endl;

// return t;
// }

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  constexpr auto CauchyGreen = StrainTags::rightCauchyGreenTensor;

  auto c0 = Eigen::Matrix3d::Identity().eval();

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
  auto bk = makeBlatzKo(ShearModulus{mu});

  std::array<double, 1> mu_og    = {mu};
  std::array<double, 1> alpha_og = {2.0};
  auto ogdenTotal                = makeOgden<1, PrincipalStretchTag::total>(mu_og, alpha_og, {Lambda}, VF3{});
  auto ogdenDevi                 = makeOgden<1, PrincipalStretchTag::deviatoric>(mu_og, alpha_og, {K}, VF3{});

  t.subTest(testMaterial<CauchyGreen>(nh, c0));
  t.subTest(testMaterial<CauchyGreen>(nh, c));

  // t.subTest(testMaterial<CauchyGreen>(bk, c0));
  t.subTest(testMaterial<CauchyGreen>(bk, c));

  // t.subTest(testMaterial<CauchyGreen>(ogdenTotal, c0));
  t.subTest(testMaterial<CauchyGreen>(ogdenTotal, c));

  // t.subTest(testMaterial<CauchyGreen>(ogdenDevi, c0));
  t.subTest(testMaterial<CauchyGreen>(ogdenDevi, c));

  // testMaterialByAD();

  return t.exit();
}