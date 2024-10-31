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
struct OgdenT
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
   * \brief Constructor for OgdenT.
   * \param mpt The Lame's parameters (first parameter and shear modulus).
   */
  explicit OgdenT(const MaterialParameters& mpt, const OgdenParameters& opt)
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

    for (auto i : parameterRange()) {
      energy += mu[i] / alpha[i] *
                (std::pow(lambdas[0], alpha[i]) + std::pow(lambdas[1], alpha[i]) + std::pow(lambdas[2], alpha[i]) - 3);
    }

    return energy;
  }

  FirstDerivative firstDerivativeImpl(const PrincipalStretches& lambdas) const {
    auto P = FirstDerivative::Zero().eval();

    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;

    // for (auto j : parameterRange())
    //   for (auto k : dimensionRange()) {
    //     auto sum = std::pow(lambdas[0], alpha[j]) + std::pow(lambdas[1], alpha[j]) + std::pow(lambdas[2], alpha[j]);
    //     // P[k] += mu[j] * std::pow(lambdas[k], alpha[j] - 1); // derivative
    //     P[k] += mu[j] * (std::pow(lambdas[k], alpha[j]) - 1.0 / 3.0 * sum ); // derivative according to
    //     Tarun MA
    //   }

    for (auto j : parameterRange()) { // iterate over Ogden terms
      auto p = (1.0 / 3.0) *
               (std::pow(lambdas[0], alpha[j]) + std::pow(lambdas[1], alpha[j]) + std::pow(lambdas[2], alpha[j]));
      for (auto k : dimensionRange()) { // iterate over principal stretches
        P[k] += mu[j] * (std::pow(lambdas[k], alpha[j]) - p);
      }
    }
    // P *= 2 / (lambdas[0] * lambdas[1] * lambdas[2]);
    return P;
  }

  SecondDerivative secondDerivativeImpl(const PrincipalStretches& lambdas) const {
    auto dS = SecondDerivative::Zero().eval();

    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;

    for (auto j : parameterRange())
      for (auto k : dimensionRange()) {
        auto p = (1.0 / 3.0) *
                 (std::pow(lambdas[0], alpha[j]) + std::pow(lambdas[1], alpha[j]) + std::pow(lambdas[2], alpha[j]));

        // dS(k, k) += mu[j] * ((alpha[j] - 1) * std::pow(lambdas[k], alpha[j] - 2) - p);
        dS(k,k) = mu[j] * alpha[j] * (alpha[j] - 1) * std::pow(lambdas[k], alpha[j] - 2) + p;
      }
    return dS;
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return OgdenT<ScalarTypeOther> The rebound Ogden material.
   */
  template <typename STO>
  auto rebind() const {
    return OgdenT<STO, nOgdenParameters>(materialParameters_);
  }

private:
  MaterialParameters materialParameters_;
  OgdenParameters ogdenParameters_;

  inline auto parameterRange() const { return Dune::Hybrid::integralRange(nOgdenParameters); }
  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(worldDimension); }
};

/**
 * \brief Alias for OgdenT with double as the default scalar type.
 */
template <int n>
using Ogden = OgdenT<double, n>;

} // namespace Ikarus
