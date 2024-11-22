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
#include <ikarus/experimental/autodiffmat.hh>

using namespace Ikarus;
using namespace Ikarus::Materials;
using namespace Ikarus::Experimental;
using Dune::TestSuite;

template <typename MAT, StrainTags strainTag>
auto checkMaterialByAutoDiffImpl(const MAT& mat, const auto C, const std::string& testName = "") {
  Dune::TestSuite t("Check storedEnergyImpl() and stressesImpl() by Automatic Differentiation" + testName);
  const auto matAD1 = AutoDiffMAT<MAT, true>(mat);
  const auto matAD2 = AutoDiffMAT<MAT, true, true>(mat);

  double energy   = mat.template storedEnergy<strainTag>(C);
  auto stress     = mat.template stresses<strainTag>(C);
  auto matTangent = mat.template tangentModuli<strainTag>(C);

  auto energy_ad1     = matAD1.template storedEnergy<strainTag>(C);
  auto stress_ad1     = matAD1.template stresses<strainTag>(C);
  auto matTangent_ad1 = matAD1.template tangentModuli<strainTag>(C);

  auto energy_ad2     = matAD2.template storedEnergy<strainTag>(C);
  auto stress_ad2     = matAD2.template stresses<strainTag>(C);
  auto matTangent_ad2 = matAD2.template tangentModuli<strainTag>(C);

  constexpr double tol = 1e-10;

  checkScalars(t, energy, static_cast<double>(energy_ad1), testLocation() + "Incorrect Energies.", tol);
  checkScalars(t, energy, static_cast<double>(energy_ad2), testLocation() + "Incorrect Energies.", tol);

  checkApproxVectors(t, stress, stress_ad1, testLocation() + "Incorrect stresses.", tol);
  checkApproxVectors(t, stress, stress_ad2, testLocation() + "Incorrect stresses.", tol);

  checkApproxMatrices(t, matTangent, matTangent_ad1, testLocation() + "Incorrect tangentModuli.", tol);
  checkApproxMatrices(t, matTangent, matTangent_ad2, testLocation() + "Incorrect tangentModuli.", tol);

  return t;
}

template <typename MAT>
auto checkMaterialByAutoDiff(const MAT& mat) {
  Dune::TestSuite t("AutoDiff Test: " + mat.name());

  constexpr StrainTags CauchyGreen = StrainTags::rightCauchyGreenTensor;

  auto C0 = Eigen::Matrix3d::Identity().eval();

  if constexpr (requires { mat.deviatoricFunction(); } and requires { mat.volumetricFunction(); }) {
    const std::string autodiffErrorMsg = "AutoDiff with duplicate principal stretches should have failed here.";
    t.checkThrow<Dune::InvalidStateException>(
        [&]() { checkMaterialByAutoDiffImpl<MAT, CauchyGreen>(mat, C0, " Undeformed State"); },
        testLocation() + autodiffErrorMsg);
  } else {
    t.subTest(checkMaterialByAutoDiffImpl<MAT, CauchyGreen>(mat, C0, " Undeformed State"));
  }

  Eigen::Matrix3d C{
      { 0.600872, -0.179083, 0},
      {-0.179083,  0.859121, 0},
      {        0,         0, 1}
  };

  t.subTest(checkMaterialByAutoDiffImpl<MAT, CauchyGreen>(mat, C, " Deformed State"));
  return t;
}