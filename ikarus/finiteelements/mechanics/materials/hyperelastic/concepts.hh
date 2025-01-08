// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file concepts.hh
 * \brief Header file including concepts for hyperelastic material models.
 * \ingroup  materials
 */

#pragma once

#include <concepts>
#include <string>

namespace Ikarus::Concepts {

/**
 * \concept DeviatoricFunction
 * \brief Concept to check if the underlying function is a deviatoric function.
 * \tparam Type of the deviatoric function.
 *
 */
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
  { dm.materialParametersImpl() } -> std::same_as<typename DF::MaterialParameters>;
  { dm.name() } -> std::convertible_to<std::string>;
};

/**
 * \concept VolumetricFunction
 * \brief Concept to check if the underlying function is a volumetric function.
 * \tparam Type of the volumetric function.
 *
 */
template <typename VF>
concept VolumetricFunction = requires(VF vf, const typename VF::JType& j) {
  typename VF::ScalarType;
  typename VF::JType;

  { vf.storedEnergyImpl(j) } -> std::same_as<typename VF::ScalarType>;
  { vf.firstDerivativeImpl(j) } -> std::same_as<typename VF::ScalarType>;
  { vf.secondDerivativeImpl(j) } -> std::same_as<typename VF::ScalarType>;
  { vf.name() } -> std::convertible_to<std::string>;
};

/**
 * \concept DeviatoricPart
 * \brief Concept to check if the underlying type corresponds to the deviatoric part of a hyperelastic material model.
 * \tparam Type of the deviatoric part of a hyperelastic material model.
 *
 */
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

/**
 * \concept VolumetricPart
 * \brief Concept to check if the underlying type corresponds to the volumetric part of a hyperelastic material model.
 * \tparam Type of the volumetric part of a hyperelastic material model.
 *
 */
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
