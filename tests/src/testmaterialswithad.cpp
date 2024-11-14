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

  t.check(isApproxSame(stress, stress_ad, prec)) << std::setprecision(16) << "Incorrect stresses." << " stress is\t"
                                                 << stress.transpose() << "\n stress_ad is\t" << stress_ad.transpose();
  t.check(isApproxSame(matTangent, matTangent_ad, prec))
      << std::setprecision(16) << "Incorrect tangentModuli derived from stresses." << " matTangent is\n"
      << matTangent << "\n matTangent_ad is\n"
      << matTangent_ad;
  t.check(isApproxSame(matTangent, matTangent_ad_e, prec))
      << std::setprecision(16) << "Incorrect tangentModuli derived from energy." << " matTangent is\n"
      << matTangent << "\n matTangent_ad_e is\n"
      << matTangent_ad_e;
  t.check(isApproxSame(matTangent_ad_e, matTangent_ad, prec))
      << std::setprecision(16) << "tangentModuli derived from energy and stresses are not equal."
      << " matTangent_ad is\n"
      << matTangent_ad << "\n matTangent_ad_e is\n"
      << matTangent_ad_e;
  ;

  return t;
}

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

  std::array<double, 3> mu_og    = {2.0 * mu / 3.0, mu / 6.0, mu / 6.0};
  std::array<double, 3> alpha_og = {1.23, 0.59, 0.18};
  auto ogdenTotal                = makeOgden<3, PrincipalStretchTag::total>(mu_og, alpha_og, {Lambda}, VF3{});
  auto ogdenDevi                 = makeOgden<3, PrincipalStretchTag::deviatoric>(mu_og, alpha_og, {K}, VF3{});
  auto mr                        = makeMooneyRivlin({mu / 2.0, mu / 2.0}, {K}, VF3{});
  auto yeoh                      = makeYeoh({2.0 * mu / 3.0, mu / 6.0, mu / 6.0}, {K}, VF3{});

  const std::string autodiffErrorMsg = "AutoDiff with duplicate principal stretches should have failed here.";

  t.subTest(testMaterial<CauchyGreen>(nh, c0));
  t.subTest(testMaterial<CauchyGreen>(nh, c));

  t.checkThrow<Dune::InvalidStateException>([&]() { testMaterial<CauchyGreen>(bk, c0); },
                                            testLocation() + autodiffErrorMsg);
  t.subTest(testMaterial<CauchyGreen>(bk, c));

  t.checkThrow<Dune::InvalidStateException>([&]() { testMaterial<CauchyGreen>(ogdenTotal, c0); },
                                            testLocation() + autodiffErrorMsg);
  t.subTest(testMaterial<CauchyGreen>(ogdenTotal, c));

  t.checkThrow<Dune::InvalidStateException>([&]() { testMaterial<CauchyGreen>(ogdenDevi, c0); },
                                            testLocation() + autodiffErrorMsg);
  t.subTest(testMaterial<CauchyGreen>(ogdenDevi, c));

  t.subTest(testMaterial<CauchyGreen>(mr, c));
  t.subTest(testMaterial<CauchyGreen>(yeoh, c));

  return t.exit();
}
