// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file concepts.hh
 * \brief Header file for material models in Ikarus finite element mechanics.
 * \ingroup  mechanics
 */

#pragma once

#include <concepts>
#include <string>

namespace Ikarus::Concepts {

template <typename DF>
concept DeviatoricFunction = requires(DF dm, const typename DF::PrincipalStretches& lambda) {
  typename DF::ScalarType;
  typename DF::PrincipalStretches;

  typename DF::FirstDerivative;
  typename DF::SecondDerivative;
  typename DF::MaterialParameters;

  { dm.storedEnergyImpl(lambda) } -> std::same_as<typename DF::ScalarType>;
  { dm.firstDerivativeImpl(lambda) } -> std::same_as<typename DF::FirstDerivative>;
  { dm.secondDerivativeImpl(lambda) } -> std::same_as<typename DF::SecondDerivative>;
};

template <typename VF>
concept VolumetricFunction = requires(VF vf, const typename VF::JType& j) {
  typename VF::ScalarType;
  typename VF::JType;

  { vf.storedEnergyImpl(j) } -> std::same_as<typename VF::ScalarType>;
  { vf.firstDerivativeImpl(j) } -> std::same_as<typename VF::ScalarType>;
  { vf.secondDerivativeImpl(j) } -> std::same_as<typename VF::ScalarType>;
};

template <typename DP>
concept DeviatoricPart = requires(DP dp, const typename DP::PrincipalStretches& lambda) {
  typename DP::DeviatoricFunction;
  typename DP::ScalarType;
  typename DP::PrincipalStretches;
  typename DP::MaterialParameters;

  typename DP::StressMatrix;
  typename DP::MaterialTensor;

  { DP::name() } -> std::convertible_to<std::string>;

  { dp.storedEnergy(lambda) } -> std::same_as<typename DP::ScalarType>;
  { dp.stresses(lambda) } -> std::same_as<typename DP::StressMatrix>;
  { dp.tangentModuli(lambda) } -> std::same_as<typename DP::MaterialTensor>;
};

template <typename VP>
concept VolumetricPart = requires(VP vp, const typename VP::JType& j) {
  typename VP::VolumetricFunction;
  typename VP::ScalarType;
  typename VP::JType;
  typename VP::MaterialParameter;

  { VP::name() } -> std::convertible_to<std::string>;

  { vp.storedEnergy(j) } -> std::same_as<typename VP::ScalarType>;
  { vp.firstDerivative(j) } -> std::same_as<typename VP::ScalarType>;
  { vp.secondDerivative(j) } -> std::same_as<typename VP::ScalarType>;
};

} // namespace Ikarus::Concepts
