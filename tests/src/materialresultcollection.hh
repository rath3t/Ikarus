// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <Eigen/Core>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/utils/makeenum.hh>

using namespace Ikarus;
using namespace Ikarus::Materials;

MAKE_ENUM(DeformationType, Undeformed, UniaxialTensile, BiaxialTensile, PureShear, Random);

template <typename ST, typename FD, typename SD>
auto initializeMaterialResults() {
  return std::make_tuple(ST{0.0}, FD::Zero().eval(), SD::Zero().eval());
}

template <DeformationType def>
auto BlatzKoResults() {
  using DEV                              = BlatzKo;
  using ST                               = typename DEV::ScalarType;
  using FD                               = typename DEV::FirstDerivative;
  using SD                               = typename DEV::SecondDerivative;
  auto [energy, stresses, tangentModuli] = initializeMaterialResults<ST, FD, SD>();

  if constexpr (def == DeformationType::UniaxialTensile) {
  } else if constexpr (def == DeformationType::BiaxialTensile) {
  } else if constexpr (def == DeformationType::PureShear) {
  } else if constexpr (def == DeformationType::UniaxialTensile) {
  } else {
  }

  return std::make_tuple(energy, stresses, tangentModuli);
}

template <DeformationType def>
auto OgdenTotalResults() {
  using DEV                              = Ogden<3, PrincipalStretchTag::total>;
  using ST                               = typename DEV::ScalarType;
  using FD                               = typename DEV::FirstDerivative;
  using SD                               = typename DEV::SecondDerivative;
  auto [energy, stresses, tangentModuli] = initializeMaterialResults<ST, FD, SD>();

  if constexpr (def == DeformationType::UniaxialTensile) {
  } else if constexpr (def == DeformationType::BiaxialTensile) {
  } else if constexpr (def == DeformationType::PureShear) {
  } else if constexpr (def == DeformationType::UniaxialTensile) {
  } else {
  }

  return std::make_tuple(energy, stresses, tangentModuli);
}

template <DeformationType def>
auto OgdenDeviatoricResults() {
  using DEV                              = Ogden<3, PrincipalStretchTag::deviatoric>;
  using ST                               = typename DEV::ScalarType;
  using FD                               = typename DEV::FirstDerivative;
  using SD                               = typename DEV::SecondDerivative;
  auto [energy, stresses, tangentModuli] = initializeMaterialResults<ST, FD, SD>();

  if constexpr (def == DeformationType::UniaxialTensile) {
  } else if constexpr (def == DeformationType::BiaxialTensile) {
  } else if constexpr (def == DeformationType::PureShear) {
  } else if constexpr (def == DeformationType::UniaxialTensile) {
  } else {
  }

  return std::make_tuple(energy, stresses, tangentModuli);
}

template <DeformationType def>
auto MooneyRivlinResults() {
  using DEV                              = InvariantBased<2>;
  using ST                               = typename DEV::ScalarType;
  using FD                               = typename DEV::FirstDerivative;
  using SD                               = typename DEV::SecondDerivative;
  auto [energy, stresses, tangentModuli] = initializeMaterialResults<ST, FD, SD>();

  if constexpr (def == DeformationType::UniaxialTensile) {
  } else if constexpr (def == DeformationType::BiaxialTensile) {
  } else if constexpr (def == DeformationType::PureShear) {
  } else if constexpr (def == DeformationType::UniaxialTensile) {
  } else {
  }

  return std::make_tuple(energy, stresses, tangentModuli);
}

template <DeformationType def>
auto YeohResults() {
  using DEV                              = InvariantBased<3>;
  using ST                               = typename DEV::ScalarType;
  using FD                               = typename DEV::FirstDerivative;
  using SD                               = typename DEV::SecondDerivative;
  auto [energy, stresses, tangentModuli] = initializeMaterialResults<ST, FD, SD>();

  if constexpr (def == DeformationType::UniaxialTensile) {
  } else if constexpr (def == DeformationType::BiaxialTensile) {
  } else if constexpr (def == DeformationType::PureShear) {
  } else if constexpr (def == DeformationType::UniaxialTensile) {
  } else {
  }

  return std::make_tuple(energy, stresses, tangentModuli);
}

template <DeformationType def>
auto ArrudaBoyceResults() {
  using DEV                              = ArrudaBoyce;
  using ST                               = typename DEV::ScalarType;
  using FD                               = typename DEV::FirstDerivative;
  using SD                               = typename DEV::SecondDerivative;
  auto [energy, stresses, tangentModuli] = initializeMaterialResults<ST, FD, SD>();

  if constexpr (def == DeformationType::UniaxialTensile) {
  } else if constexpr (def == DeformationType::BiaxialTensile) {
  } else if constexpr (def == DeformationType::PureShear) {
  } else if constexpr (def == DeformationType::UniaxialTensile) {
  } else {
  }

  return std::make_tuple(energy, stresses, tangentModuli);
}

template <DeformationType def>
auto GentResults() {
  using DEV                              = Gent;
  using ST                               = typename DEV::ScalarType;
  using FD                               = typename DEV::FirstDerivative;
  using SD                               = typename DEV::SecondDerivative;
  auto [energy, stresses, tangentModuli] = initializeMaterialResults<ST, FD, SD>();

  if constexpr (def == DeformationType::UniaxialTensile) {
  } else if constexpr (def == DeformationType::BiaxialTensile) {
  } else if constexpr (def == DeformationType::PureShear) {
  } else if constexpr (def == DeformationType::UniaxialTensile) {
  } else {
  }

  return std::make_tuple(energy, stresses, tangentModuli);
}