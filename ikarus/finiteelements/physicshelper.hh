// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file physicshelper.hh
 * \brief Material property functions and conversion utilities.
 * \details This file contains material property functions and conversion utilities related to linear elasticity.
 */

#pragma once
#include <dune/common/float_cmp.hh>

#include <Eigen/Core>
namespace Ikarus {

///< Structure representing Young's modulus and Poisson's ratio. \brief see
///< https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
struct YoungsModulusAndPoissonsRatio
{
  double emodul;
  double nu;
};

///< Structure representing Young's modulus and shear modulus.
struct YoungsModulusAndShearModulus
{
  double emodul;
  double mu;
};

///< Structure representing Young's modulus and bulk modulus.
struct YoungsModulusAndBulkModulus
{
  double emodul;
  double K;
};

///< Structure representing Young's modulus and Lame's first parameter.
struct YoungsModulusAndLamesFirstParameter
{
  double emodul;
  double lambda;
};

///< Structure representing bulk modulus and Lame's first parameter.
struct BulkModulusAndLamesFirstParameter
{
  double K;
  double lambda;
};

///< Structure representing Lame's first parameter and shear modulus.
struct LamesFirstParameterAndShearModulus
{
  double lambda;
  double mu;
};

/**
 * \brief Concept for checking if a type is a valid material parameter tuple.
 *
 * \tparam MP Type to check.
 */
template <typename MP>
concept MPTuple =
    std::is_same_v<MP, YoungsModulusAndPoissonsRatio> or std::is_same_v<MP, YoungsModulusAndBulkModulus> or
    std::is_same_v<MP, YoungsModulusAndLamesFirstParameter> or std::is_same_v<MP, BulkModulusAndLamesFirstParameter> or
    std::is_same_v<MP, LamesFirstParameterAndShearModulus> or std::is_same_v<MP, YoungsModulusAndShearModulus>;

/**
 * \brief Conversion utility for Lame's constants.
 *
 * \tparam ValuePair Type of the value pair to convert.
 */
template <typename ValuePair>
struct ConvertLameConstants
{
  constexpr double toLamesFirstParameter()
  requires(!std::is_same_v<ValuePair, YoungsModulusAndLamesFirstParameter> and
           !std::is_same_v<ValuePair, BulkModulusAndLamesFirstParameter> and
           !std::is_same_v<ValuePair, LamesFirstParameterAndShearModulus>)
  {
    if constexpr (std::is_same_v<ValuePair, YoungsModulusAndPoissonsRatio>) {
      const auto& E  = vp_.emodul;
      const auto& nu = vp_.nu;
      return Dune::FloatCmp::eq(nu, 0.5) ? std::numeric_limits<double>::infinity()
                                         : E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    } else if constexpr (std::is_same_v<ValuePair, YoungsModulusAndShearModulus>) {
      const auto& E  = vp_.emodul;
      const auto& mu = vp_.mu;
      return mu * (E - 2.0 * mu) / (3.0 * mu - E);
    } else if constexpr (std::is_same_v<ValuePair, YoungsModulusAndBulkModulus>) {
      const auto& E = vp_.emodul;
      const auto& K = vp_.K;
      return 3.0 * K * (3.0 * K - E) / (9.0 * K - E);
    } else
      assert(false && "Your LameParameter request is not implemented");
  }

  constexpr double toBulkModulus()
  requires(!std::is_same_v<ValuePair, YoungsModulusAndBulkModulus> and
           !std::is_same_v<ValuePair, BulkModulusAndLamesFirstParameter>)
  {
    if constexpr (std::is_same_v<ValuePair, YoungsModulusAndPoissonsRatio>) {
      const auto& E  = vp_.emodul;
      const auto& nu = vp_.nu;
      return E / (3.0 * (1.0 - 2.0 * nu));
    } else if constexpr (std::is_same_v<ValuePair, YoungsModulusAndShearModulus>) {
      const auto& E  = vp_.emodul;
      const auto& mu = vp_.mu;
      return E * mu / (3.0 * (3.0 * mu - E));
    } else if constexpr (std::is_same_v<ValuePair, YoungsModulusAndLamesFirstParameter>) {
      const auto& E      = vp_.emodul;
      const auto& lambda = vp_.lambda;
      return (E + 3.0 * lambda + calcR(vp_)) / 6.0;
    } else if constexpr (std::is_same_v<ValuePair, LamesFirstParameterAndShearModulus>) {
      const auto& lambda = vp_.lambda;
      const auto& mu     = vp_.mu;
      return lambda + 2.0 * mu / 3.0;
    } else
      assert(false && "Your LameParameter request is not implemented");
  }

  constexpr double toShearModulus()
  requires(!std::is_same_v<ValuePair, YoungsModulusAndShearModulus> and
           !std::is_same_v<ValuePair, LamesFirstParameterAndShearModulus>)
  {
    if constexpr (std::is_same_v<ValuePair, YoungsModulusAndPoissonsRatio>) {
      const auto& E  = vp_.emodul;
      const auto& nu = vp_.nu;
      return E / (2.0 * (1.0 + nu));
    } else if constexpr (std::is_same_v<ValuePair, YoungsModulusAndBulkModulus>) {
      const auto& E = vp_.emodul;
      const auto& K = vp_.K;
      return 3.0 * K * E / (9.0 * K - E);
    } else if constexpr (std::is_same_v<ValuePair, YoungsModulusAndLamesFirstParameter>) {
      const auto& E      = vp_.emodul;
      const auto& lambda = vp_.lambda;
      return (E - 3.0 * lambda + calcR(vp_)) / 4.0;
    } else if constexpr (std::is_same_v<ValuePair, BulkModulusAndLamesFirstParameter>) {
      const auto& K      = vp_.K;
      const auto& lambda = vp_.lambda;
      return 3.0 * (K - lambda) / 2.0;
    } else
      assert(false && "Your LameParameter request is not implemented");
  }

  constexpr double toPWaveModulus() {
    if constexpr (std::is_same_v<ValuePair, YoungsModulusAndPoissonsRatio>) {
      const auto& E  = vp_.emodul;
      const auto& nu = vp_.nu;
      return E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    } else if constexpr (std::is_same_v<ValuePair, YoungsModulusAndShearModulus>) {
      const auto& E  = vp_.emodul;
      const auto& mu = vp_.mu;
      return mu * (4.0 * mu - E) / (3.0 * mu - E);
    } else if constexpr (std::is_same_v<ValuePair, YoungsModulusAndBulkModulus>) {
      const auto& E = vp_.emodul;
      const auto& K = vp_.K;
      return 3.0 * K * (3.0 * K + E) / (9.0 * K - E);
    } else if constexpr (std::is_same_v<ValuePair, YoungsModulusAndLamesFirstParameter>) {
      const auto& E      = vp_.emodul;
      const auto& lambda = vp_.lambda;
      return (E - lambda + calcR(vp_)) / 2.0;
    } else if constexpr (std::is_same_v<ValuePair, BulkModulusAndLamesFirstParameter>) {
      const auto& K      = vp_.K;
      const auto& lambda = vp_.lambda;
      return 3.0 * K - 2.0 * lambda;
    } else if constexpr (std::is_same_v<ValuePair, LamesFirstParameterAndShearModulus>) {
      const auto& lambda = vp_.lambda;
      const auto& mu     = vp_.mu;
      return lambda + 2.0 * mu;
    } else
      assert(false && "Your LameParameter request is not implemented");
  }

  constexpr double toPoissonsRatio()
  requires(!std::is_same_v<ValuePair, YoungsModulusAndPoissonsRatio>)
  {
    if constexpr (std::is_same_v<ValuePair, YoungsModulusAndShearModulus>) {
      const auto& E  = vp_.emodul;
      const auto& mu = vp_.mu;
      return E / (2.0 * mu) - 1.0;
    } else if constexpr (std::is_same_v<ValuePair, YoungsModulusAndBulkModulus>) {
      const auto& E = vp_.emodul;
      const auto& K = vp_.K;
      return (3.0 * K - E) / (6.0 * K);
    } else if constexpr (std::is_same_v<ValuePair, YoungsModulusAndLamesFirstParameter>) {
      const auto& E      = vp_.emodul;
      const auto& lambda = vp_.lambda;
      return 2.0 * lambda / (E + lambda + calcR(vp_));
    } else if constexpr (std::is_same_v<ValuePair, BulkModulusAndLamesFirstParameter>) {
      const auto& K      = vp_.K;
      const auto& lambda = vp_.lambda;
      return lambda / (3 * K - lambda);
    } else if constexpr (std::is_same_v<ValuePair, LamesFirstParameterAndShearModulus>) {
      const auto& lambda = vp_.lambda;
      const auto& mu     = vp_.mu;
      return lambda / (2.0 * (lambda + mu));
    } else
      assert(false && "Your LameParameter request is not implemented");
  }

  constexpr double toYoungsModulus()
  requires(!std::is_same_v<ValuePair, YoungsModulusAndPoissonsRatio> and
           !std::is_same_v<ValuePair, YoungsModulusAndShearModulus> and
           !std::is_same_v<ValuePair, YoungsModulusAndBulkModulus> and
           !std::is_same_v<ValuePair, YoungsModulusAndLamesFirstParameter>)
  {
    if constexpr (std::is_same_v<ValuePair, BulkModulusAndLamesFirstParameter>) {
      return 9.0 * vp_.K * (vp_.K - vp_.lambda) / (3.0 * vp_.K - vp_.lambda);
    } else if constexpr (std::is_same_v<ValuePair, LamesFirstParameterAndShearModulus>) {
      const auto& lambda = vp_.lambda;
      const auto& mu     = vp_.mu;
      return mu * (3.0 * lambda + 2.0 * mu) / (lambda + mu);
    } else
      assert(false && "Your LameParameter request is not implemented");
  }

private:
  friend ConvertLameConstants<YoungsModulusAndPoissonsRatio> convertLameConstants(
      const YoungsModulusAndPoissonsRatio& valuePair);
  friend ConvertLameConstants<YoungsModulusAndShearModulus> convertLameConstants(
      const YoungsModulusAndShearModulus& valuePair);

  friend ConvertLameConstants<YoungsModulusAndBulkModulus> convertLameConstants(
      const YoungsModulusAndBulkModulus& valuePair);

  friend ConvertLameConstants<LamesFirstParameterAndShearModulus> convertLameConstants(
      const LamesFirstParameterAndShearModulus& valuePair);

  friend ConvertLameConstants<BulkModulusAndLamesFirstParameter> convertLameConstants(
      const BulkModulusAndLamesFirstParameter& valuePair);

  friend ConvertLameConstants<YoungsModulusAndLamesFirstParameter> convertLameConstants(
      const YoungsModulusAndLamesFirstParameter& valuePair);

  ConvertLameConstants(ValuePair&& valuePair)
      : vp_(valuePair) {}
  ConvertLameConstants(const ValuePair& valuePair)
      : vp_(valuePair) {}

  double calcR(const YoungsModulusAndLamesFirstParameter& valuePair) {
    const auto& E      = valuePair.emodul;
    const auto& lambda = valuePair.lambda;
    return std::sqrt(E * E + 9 * lambda * lambda + 2 * E * lambda);
  }
  ValuePair vp_;
};
inline ConvertLameConstants<YoungsModulusAndPoissonsRatio> convertLameConstants(
    const YoungsModulusAndPoissonsRatio& valuePair) {
  return {valuePair};
}
inline ConvertLameConstants<YoungsModulusAndShearModulus> convertLameConstants(
    const YoungsModulusAndShearModulus& valuePair) {
  return {valuePair};
}
inline ConvertLameConstants<YoungsModulusAndBulkModulus> convertLameConstants(
    const YoungsModulusAndBulkModulus& valuePair) {
  return {valuePair};
}
inline ConvertLameConstants<LamesFirstParameterAndShearModulus> convertLameConstants(
    const LamesFirstParameterAndShearModulus& valuePair) {
  return {valuePair};
}
inline ConvertLameConstants<BulkModulusAndLamesFirstParameter> convertLameConstants(
    const BulkModulusAndLamesFirstParameter& valuePair) {
  return {valuePair};
}
inline ConvertLameConstants<YoungsModulusAndLamesFirstParameter> convertLameConstants(
    const YoungsModulusAndLamesFirstParameter& valuePair) {
  return {valuePair};
}

/**
 * \brief Converts Young's modulus and Poisson's ratio to Lame's first parameter and shear modulus.
 *
 * \param matParameter Young's modulus and Poisson's ratio.
 * \return Lame's first parameter and shear modulus.
 */
inline auto toLamesFirstParameterAndShearModulus(const YoungsModulusAndPoissonsRatio& matParameter) {
  auto lambda = convertLameConstants(matParameter).toLamesFirstParameter();
  auto mu     = convertLameConstants(matParameter).toShearModulus();

  return LamesFirstParameterAndShearModulus{.lambda = lambda, .mu = mu};
}

/**
 * \brief Converts Lame's first parameter and shear modulus to Young's modulus and Poisson's ratio.
 *
 * \param matParameter Lame's first parameter and shear modulus.
 * \return Young's modulus and Poisson's ratio.
 */
inline auto toYoungsModulusAndPoissonsRatio(const LamesFirstParameterAndShearModulus& matParameter) {
  auto emod = convertLameConstants(matParameter).toYoungsModulus();
  auto nu   = convertLameConstants(matParameter).toPoissonsRatio();

  return YoungsModulusAndPoissonsRatio{.emodul = emod, .nu = nu};
}

} // namespace Ikarus
