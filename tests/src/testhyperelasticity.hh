// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "materialresultcollection.hh"
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

auto recoverNeoHookeTest() {
  using StrainTags::rightCauchyGreenTensor;

  TestSuite t("Recover NeoHooke material");

  auto matPar = testMatPar();
  auto mu     = convertLameConstants(matPar).toShearModulus();
  auto Lambda = convertLameConstants(matPar).toLamesFirstParameter();
  auto c      = testMatrix();
  auto K      = convertLameConstants(matPar).toBulkModulus();

  std::array<double, 1> mu_og    = {mu};
  std::array<double, 1> alpha_og = {2.0};

  auto nhFromogdenTotal     = makeOgden<1, PrincipalStretchTag::total>(mu_og, alpha_og);
  auto nhFromogdenDevi      = makeOgden<1, PrincipalStretchTag::deviatoric>(mu_og, alpha_og);
  auto nhFromInvariantBased = makeInvariantBased<1>({mu / 2.0}, {1}, {0});
  auto nh                   = NeoHooke(toLamesFirstParameterAndShearModulus(matPar));

  constexpr double tol = 1e-14;

  auto checkNHRecovery = [&]<typename MAT1, typename MAT2>(const MAT1& mat1, const MAT2& mat2) {
    auto energy_mat1 = mat1.template storedEnergy<rightCauchyGreenTensor>(c);
    auto stress_mat1 = mat1.template stresses<rightCauchyGreenTensor>(c);
    auto moduli_mat1 = mat1.template tangentModuli<rightCauchyGreenTensor>(c);

    auto energy_mat2 = mat2.template storedEnergy<rightCauchyGreenTensor>(c);
    auto stress_mat2 = mat2.template stresses<rightCauchyGreenTensor>(c);
    auto moduli_mat2 = mat2.template tangentModuli<rightCauchyGreenTensor>(c);

    const std::string matName = mat1.name() + " and " + mat2.name();

    checkScalars(t, energy_mat2, energy_mat1, testLocation() + matName + "Incorrect Energy.", tol);
    checkApproxVectors(t, stress_mat2, stress_mat1, testLocation() + matName + "Incorrect Stresses.", tol);
    checkApproxMatrices(t, moduli_mat2, moduli_mat1, testLocation() + matName + "Incorrect tangentModuli.", tol);
  };

  checkNHRecovery(nh, nhFromogdenTotal);
  checkNHRecovery(nhFromInvariantBased, nhFromogdenDevi);

  return t;
}

struct Deformations
{
  static constexpr int dim = 3;
  using MatrixType         = Eigen::Matrix<double, dim, dim>;
  Deformations()           = default;

  template <DeformationType DT>
  MatrixType rightCauchyGreen(double lambda) const {
    if constexpr (DT == DeformationType::Undeformed)
      return undeformed(lambda);
    else if constexpr (DT == DeformationType::UniaxialTensile)
      return uniaxialTensile(lambda);
    else if constexpr (DT == DeformationType::BiaxialTensile)
      return biaxialTensile(lambda);
    else if constexpr (DT == DeformationType::PureShear)
      return pureShear(lambda);
    else
      return random(lambda);
  }

private:
  /// Convert deformation gradient to right Cauchy-Green tensor
  MatrixType toC(const MatrixType& F) const { return (F.transpose() * F).eval(); }

  MatrixType undeformed([[maybe_unused]] double lambda_) const {
    auto F = MatrixType::Identity().eval();
    return toC(F);
  }

  MatrixType uniaxialTensile(double lambda_) const {
    auto F = MatrixType::Zero().eval();
    F.diagonal() << lambda_, 1.0 / sqrt(lambda_), 1.0 / sqrt(lambda_);
    return toC(F);
  }

  MatrixType biaxialTensile(double lambda_) const {
    auto F = MatrixType::Zero().eval();
    F.diagonal() << lambda_, lambda_, 1.0 / (lambda_ * lambda_);
    return toC(F);
  }

  MatrixType pureShear(double lambda_) const {
    auto F = MatrixType::Zero().eval();
    F.diagonal() << lambda_, 1.0, 1.0 / lambda_;
    return toC(F);
  }

  MatrixType random([[maybe_unused]] double lambda_) const { return testMatrix(); }
};

template <typename DEV, DeformationType def>
requires(Concepts::DeviatoricFunction<DEV>)
auto materialResults() {
  if constexpr (std::is_same_v<DEV, BlatzKoT<double>>) {
    return BlatzKoResults<def>();
  } else if constexpr (std::is_same_v<DEV, OgdenT<double, 3, PrincipalStretchTag::total>>) {
    return OgdenTotalResults<def>();
  } else if constexpr (std::is_same_v<DEV, OgdenT<double, 3, PrincipalStretchTag::deviatoric>>) {
    return OgdenDeviatoricResults<def>();
  } else if constexpr (std::is_same_v<DEV, InvariantBasedT<double, 2>>) {
    return MooneyRivlinResults<def>();
  } else if constexpr (std::is_same_v<DEV, InvariantBasedT<double, 3>>) {
    return YeohResults<def>();
  } else if constexpr (std::is_same_v<DEV, ArrudaBoyceT<double>>) {
    return ArrudaBoyceResults<def>();
  } else if constexpr (std::is_same_v<DEV, GentT<double>>) {
    return GentResults<def>();
  } else
    static_assert(Dune::AlwaysFalse<DEV>::value, "The requested deviatoric function is not implemented.");
}

template <DeformationType def, typename DEV, StrainTags tag = StrainTags::rightCauchyGreenTensor>
requires(Concepts::DeviatoricFunction<DEV>)
auto testMaterialResult(const DEV& dev) {
  Dune::TestSuite t("Test Deviatoric Function Results for the material model: " + dev.name() +
                    " with deformation type as " + toString(def));
  auto [energy_ex, stresses_ex, tangentModuli_ex] = materialResults<DEV, def>();
  auto deformation                                = Deformations{};
  constexpr double lambda                         = 1.37;
  auto C = Materials::Impl::maybeFromVoigt(deformation.rightCauchyGreen<def>(lambda));
  Eigen::SelfAdjointEigenSolver<decltype(C)> eigensolver{};
  eigensolver.compute(C, Eigen::EigenvaluesOnly);
  auto& principalStretches = eigensolver.eigenvalues().array().sqrt().eval();
  auto W                   = dev.storedEnergyImpl(principalStretches);
  auto dWdLambda           = dev.firstDerivativeImpl(principalStretches);
  auto ddWdLambda          = dev.secondDerivativeImpl(principalStretches);
  constexpr double tol     = 1e-14;

  checkScalars(t, W, energy_ex, testLocation() + dev.name() + ": Incorrect Energies", tol);
  checkApproxVectors(t, dWdLambda, stresses_ex, testLocation() + dev.name() + ": Incorrect stresses", tol);
  checkApproxMatrices(t, ddWdLambda, tangentModuli_ex, testLocation() + dev.name() + ": Incorrect tangentModuli", tol);

  return t;
}

auto testMaterialResults() {
  Dune::TestSuite t("Test Deviatoric Function Results ");
  auto matPar     = testMatPar();
  auto mu         = convertLameConstants(matPar).toShearModulus();
  auto K          = convertLameConstants(matPar).toBulkModulus();
  auto Lambda     = convertLameConstants(matPar).toLamesFirstParameter();
  auto lameMatPar = toLamesFirstParameterAndShearModulus(matPar);

  std::array<double, 3> mu_og    = {2.0 * mu / 3.0, mu / 6.0, mu / 6.0};
  std::array<double, 3> alpha_og = {1.23, 0.59, 0.18};

  auto bk = BlatzKo(ShearModulus{mu});

  auto ogdenTotal = Ogden<3, PrincipalStretchTag::total>(mu_og, alpha_og);
  auto ogdenDevi  = Ogden<3, PrincipalStretchTag::deviatoric>(mu_og, alpha_og);
  auto mr         = InvariantBased<2>({1, 0}, {0, 1}, {mu / 2.0, mu / 2.0});
  auto yeoh       = InvariantBased<3>({1, 2, 3}, {0, 0, 0}, {2.0 * mu / 3.0, mu / 6.0, mu / 6.0});
  auto ab         = ArrudaBoyce({mu, 0.85});
  auto gent       = Gent({mu, 0.85});

  auto checkForDeformationType = [&]<DeformationType def>() {
    t.subTest(testMaterialResult<def>(bk));
    t.subTest(testMaterialResult<def>(ogdenTotal));
    t.subTest(testMaterialResult<def>(ogdenDevi));
    t.subTest(testMaterialResult<def>(mr));
    t.subTest(testMaterialResult<def>(yeoh));
    t.subTest(testMaterialResult<def>(ab));
    t.subTest(testMaterialResult<def>(gent));
  };

  checkForDeformationType.operator()<DeformationType::Undeformed>();
  checkForDeformationType.operator()<DeformationType::UniaxialTensile>();
  checkForDeformationType.operator()<DeformationType::BiaxialTensile>();
  checkForDeformationType.operator()<DeformationType::PureShear>();
  checkForDeformationType.operator()<DeformationType::Random>();

  return t;
}
