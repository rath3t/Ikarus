// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file concepts.hh
 * \brief Header file for material models in Ikarus finite element mechanics.
 * \ingroup  mechanics
 */

#pragma once

// #include <ikarus/finiteelements/mechanics/materials/hyperelasticity/blatzko.hh>
// #include <ikarus/finiteelements/mechanics/materials/hyperelasticity/deviatoric.hh>
// #include <ikarus/finiteelements/mechanics/materials/hyperelasticity/hyperelastic.hh>
// #include <ikarus/finiteelements/mechanics/materials/hyperelasticity/modogden.hh>
// #include <ikarus/finiteelements/mechanics/materials/hyperelasticity/volumetric.hh>

#include <concepts>
#include <string>

namespace Ikarus::Concepts {

template <typename DP>
concept DeviatoricPartConcept = requires(DP dp, const typename DP::PrincipalStretches& lambda) {
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

  // rebind
};

template <typename VP>
concept VolumetricPartConcept = requires(VP vp, const typename VP::JType& j) {
  typename VP::VolumetricFunction;
  typename VP::ScalarType;
  typename VP::JType;

  { VP::name() } -> std::convertible_to<std::string>;

  { vp.storedEnergy(j) } -> std::same_as<typename VP::ScalarType>;
  { vp.firstDerivative(j) } -> std::same_as<typename VP::ScalarType>;
  { vp.secondDerivative(j) } -> std::same_as<typename VP::ScalarType>;
};

} // namespace Ikarus::Concepts
