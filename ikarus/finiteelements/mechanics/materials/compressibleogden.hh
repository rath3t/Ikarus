// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file CompressibleOgden.hh
 * \brief Implementation of the CompressibleOgden material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

template <typename ST, int n>
struct CompressibleOgdenT : public Material<CompressibleOgdenT<ST, n>>
{
  using ScalarType         = ST;
  using PrincipalStretches = Eigen::Vector<ScalarType, 3>;

  static constexpr int nOgdenParameters = n;
  static constexpr int worldDimension   = 3;

  using FirstDerivative  = Eigen::Vector<ScalarType, worldDimension>;
  using SecondDerivative = Eigen::Matrix<ScalarType, worldDimension, worldDimension>;

  using MaterialParameters = std::array<double, nOgdenParameters>;
  using OgdenParameters    = std::array<double, nOgdenParameters>;

  [[nodiscard]] constexpr static std::string nameImpl() noexcept { return "CompressibleOgden"; }

  /**
   * \brief Constructor for CompressibleOgdenT.
   * \param mpt The Lame's parameters (first parameter and shear modulus).
   */
  explicit CompressibleOgdenT(const MaterialParameters& mpt, const OgdenParameters& opt)
      : materialParameters_{mpt},
        ogdenParameters_{opt} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameters_; }

  /**
   * \brief Computes the stored energy in the Neo-Hookean material model.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return ScalarType The stored energy.
   */
  ScalarType storedEnergyImpl(const PrincipalStretches& lambdas) const {
    ScalarType energy{};
    auto J = lambdas[0] * lambdas[1] * lambdas[2];

    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;

    for (auto i : parameterRange()) {
      energy +=
          mu[i] / alpha[i] *
              (pow(lambdas[0], alpha[i]) + pow(lambdas[1], alpha[i]) + pow(lambdas[2], alpha[i]) - 3) -
          mu[i] * log(J);
    }

    return energy;
  }

  /**
   * \brief Computes the stresses in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return StressMatrix The stresses.
   */
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

  /**
   * \brief Computes the tangent moduli in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> The tangent moduli.
   */

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
   * \return CompressibleOgdenT<ScalarTypeOther> The rebound CompressibleOgden material.
   */
  template <typename STO>
  auto rebind() const {
    return CompressibleOgdenT<STO, nOgdenParameters>(materialParameters_, ogdenParameters_);
  }

private:
  MaterialParameters materialParameters_;
  OgdenParameters ogdenParameters_;

  inline auto parameterRange() const { return Dune::Hybrid::integralRange(nOgdenParameters); }
  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(worldDimension); }

};

/**
 * \brief Alias for CompressibleOgdenT with double as the default scalar type.
 */
template <int n>
using CompressibleOgden = CompressibleOgdenT<double, n>;

} // namespace Ikarus
