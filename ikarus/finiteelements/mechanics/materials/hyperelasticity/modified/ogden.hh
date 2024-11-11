// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file CompressibleOgden.hh
 * \brief Implementation of the CompressibleOgden material model.
 * \ingroup  materials
 */

#pragma once

#include <string>

#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

template <typename ST, int n>
struct ModOgdenT
{
  using ScalarType         = ST;
  using PrincipalStretches = Eigen::Vector<ScalarType, 3>;

  static constexpr int nOgdenParameters = n;
  static constexpr int dim              = 3;

  using FirstDerivative  = Eigen::Vector<ScalarType, dim>;
  using SecondDerivative = Eigen::Matrix<ScalarType, dim, dim>;

  using MaterialParameters = std::array<double, nOgdenParameters>;
  using OgdenParameters    = std::array<double, nOgdenParameters>;

  [[nodiscard]] constexpr static std::string name() noexcept {
    return "Modified Ogden (n=" + std::to_string(nOgdenParameters) + ")";
  }

  /**
   * \brief Constructor for ModifiedOgdenT.
   * \param mpt The Lame's parameters (first parameter and shear modulus).
   */
  explicit ModOgdenT(const MaterialParameters& mpt, const OgdenParameters& opt)
      : materialParameters_{mpt},
        ogdenParameters_{opt} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameters_; }

  ScalarType storedEnergyImpl(const PrincipalStretches& lambdas) const {
    ScalarType energy{};
    auto J = lambdas[0] * lambdas[1] * lambdas[2];

    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;

    for (auto i : parameterRange()) {
      energy +=
          mu[i] / alpha[i] * (pow(lambdas[0], alpha[i]) + pow(lambdas[1], alpha[i]) + pow(lambdas[2], alpha[i]) - 3) -
          mu[i] * log(J);
    }

    return energy;
  }

  FirstDerivative firstDerivativeImpl(const PrincipalStretches& lambdas) const {
    auto P = FirstDerivative::Zero().eval();

    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;

    for (auto j : parameterRange()) {
      for (auto k : dimensionRange()) {
        P[k] += (mu[j] * (pow(lambdas[k], alpha[j]) - 1)) / lambdas[k];
      }
    }
    return P;
  }

  SecondDerivative secondDerivativeImpl(const PrincipalStretches& lambdas) const {
    auto dS = SecondDerivative::Zero().eval();

    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;

    for (auto j : parameterRange())
      for (auto k : dimensionRange()) {
        dS(k, k) += -2 * (mu[j] * (pow(lambdas[k], alpha[j]) - 1)) / pow(lambdas[k], 3) +
                    (mu[j] * pow(lambdas[k], alpha[j]) * alpha[j] / lambdas[k]) / pow(lambdas[k], 2);
      }

    return dS;
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return ModifiedOgdenT<ScalarTypeOther> The rebound CompressibleOgden material.
   */
  template <typename STO>
  auto rebind() const {
    return ModOgdenT<STO, nOgdenParameters>(materialParameters_, ogdenParameters_);
  }

private:
  MaterialParameters materialParameters_;
  OgdenParameters ogdenParameters_;

  inline auto parameterRange() const { return Dune::Hybrid::integralRange(nOgdenParameters); }
  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(dim); }
};

/**
 * \brief Alias for ModifiedOgdenT with double as the default scalar type.
 */
template <int n>
using ModOgden = ModOgdenT<double, n>;

} // namespace Ikarus
