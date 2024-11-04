// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file Ogden.hh
 * \brief Implementation of the Ogden material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

template <typename ST, int n>
struct IncompressibleOgdenT
{
  using ScalarType         = ST;
  using PrincipalStretches = Eigen::Vector<ScalarType, 3>;

  static constexpr int nOgdenParameters = n;
  static constexpr int worldDimension   = 3;

  using FirstDerivative  = Eigen::Vector<ScalarType, worldDimension>;
  using SecondDerivative = Eigen::Matrix<ScalarType, worldDimension, worldDimension>;

  using MaterialParameters = std::array<double, nOgdenParameters>;
  using OgdenParameters    = std::array<double, nOgdenParameters>;

  [[nodiscard]] constexpr static std::string nameImpl() noexcept { return "Ogden"; }

  /**
   * \brief Constructor for IncompressibleOgdenT.
   * \param mpt The Lame's parameters (first parameter and shear modulus).
   */
  explicit IncompressibleOgdenT(const MaterialParameters& mpt, const OgdenParameters& opt)
      : materialParameters_{mpt},
        ogdenParameters_{opt} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameters_; }

  OgdenParameters ogdenParameters() const { return materialParameters_; }

  ScalarType storedEnergyImpl(const PrincipalStretches& lambdas) const {
    ScalarType energy{};

    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;

    for (auto i : parameterRange())
      energy +=
          mu[i] / alpha[i] * (pow(lambdas[0], alpha[i]) + pow(lambdas[1], alpha[i]) + pow(lambdas[2], alpha[i]) - 3);

    return energy;
  }

  FirstDerivative firstDerivativeImpl(const PrincipalStretches& lambdaBar) const {
    auto P = FirstDerivative::Zero().eval();

    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;
    for (auto j : parameterRange()) {
      auto p = (1.0 / 3.0) * (pow(lambdaBar[0], alpha[j]) + pow(lambdaBar[1], alpha[j]) + pow(lambdaBar[2], alpha[j]));

      for (auto k : dimensionRange()) {
        P[k] += mu[j] * (pow(lambdaBar[k], alpha[j] - 1) - p);
      }
    }

    return P;
  }

  SecondDerivative secondDerivativeImpl(const PrincipalStretches& lambdas) const {
    auto dS = SecondDerivative::Zero().eval();

    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;

    for (auto j : parameterRange())
      for (auto k : dimensionRange())
        dS(k, k) = mu[j] * alpha[j] * (alpha[j] - 1) * pow(lambdas[k], alpha[j] - 2);

    return dS;
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return IncompressibleOgdenT<ScalarTypeOther> The rebound Ogden material.
   */
  template <typename STO>
  auto rebind() const {
    return IncompressibleOgdenT<STO, nOgdenParameters>(materialParameters_, ogdenParameters_);
  }

private:
  MaterialParameters materialParameters_;
  OgdenParameters ogdenParameters_;

  inline auto parameterRange() const { return Dune::Hybrid::integralRange(nOgdenParameters); }
  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(worldDimension); }
};

/**
 * \brief Alias for IncompressibleOgdenT with double as the default scalar type.
 */
template <int n>
using IncompressibleOgden = IncompressibleOgdenT<double, n>;

} // namespace Ikarus
