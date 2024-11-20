// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file ogden.hh
 * \brief Implementation of the Ogden material model.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Implementation of the Ogden material model.
 *
 *\tparam ST The underlying scalar type.
 * \tparam n number of ogden parameters
 * \tparam tag type of principal stretch quantity, either total stretches or deviatoric stretches
 * \ingroup materials
 */
template <typename ST, int n, PrincipalStretchTag tag>
struct OgdenT
{
  using ScalarType         = ST;
  using PrincipalStretches = Eigen::Vector<ScalarType, 3>;

  static constexpr PrincipalStretchTag stretchTag = tag;
  static constexpr int nOgdenParameters           = n;
  static constexpr int dim                        = 3;
  static constexpr bool usesDeviatoricStretches   = stretchTag == PrincipalStretchTag::deviatoric;

  using FirstDerivative  = Eigen::Vector<ScalarType, dim>;
  using SecondDerivative = Eigen::Matrix<ScalarType, dim, dim>;

  using MaterialParameters = std::array<double, nOgdenParameters>;
  using OgdenParameters    = std::array<double, nOgdenParameters>;

  [[nodiscard]] constexpr static std::string name() noexcept {
    return "Ogden (n = " + std::to_string(nOgdenParameters) + ", stretch type = " + toString(tag) + ")";
  }

  /**
   * \brief Constructor for OgdenT
   *
   * \param mpt material parameters (array of mu values)
   * \param opt ogden parameters (array of alpha values)
   */
  explicit OgdenT(const MaterialParameters& mpt, const OgdenParameters& opt)
      : materialParameters_{mpt},
        ogdenParameters_{opt} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  const MaterialParameters& materialParametersImpl() const { return materialParameters_; }

  /**
   * \brief Returns the ogden parameters stored in the material
   */
  const OgdenParameters& ogdenParameters() const { return ogdenParameters_; }

  /**
   * \brief Computes the stored energy in the Ogden material model.
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  ScalarType storedEnergyImpl(const PrincipalStretches& lambda) const {
    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;

    ScalarType energy{};

    if constexpr (usesDeviatoricStretches) {
      auto lambdaBar = Impl::deviatoricStretches(lambda);

      for (auto i : parameterRange())
        energy += mu[i] / alpha[i] *
                  (pow(lambdaBar[0], alpha[i]) + pow(lambdaBar[1], alpha[i]) + pow(lambdaBar[2], alpha[i]) - 3);

    } else {
      auto J = lambda[0] * lambda[1] * lambda[2];

      for (auto i : parameterRange()) {
        energy +=
            mu[i] / alpha[i] * (pow(lambda[0], alpha[i]) + pow(lambda[1], alpha[i]) + pow(lambda[2], alpha[i]) - 3) -
            mu[i] * log(J);
      }
    }
    return energy;
  }

  /**
   * \brief Computes the first derivative of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  FirstDerivative firstDerivativeImpl(const PrincipalStretches& lambda) const {
    auto& mu       = materialParameters_;
    auto& alpha    = ogdenParameters_;
    auto dWdLambda = FirstDerivative::Zero().eval();

    if constexpr (usesDeviatoricStretches) {
      auto lambdaBar = Impl::deviatoricStretches(lambda);

      auto dWdLambdaBar = FirstDerivative::Zero().eval();
      for (auto j : parameterRange())
        for (auto k : dimensionRange())
          dWdLambdaBar[k] += mu[j] * (pow(lambdaBar[k], alpha[j] - 1));

      ScalarType sumLambdaBar{0.0};
      for (auto b : dimensionRange())
        sumLambdaBar += lambdaBar[b] * dWdLambdaBar[b];

      for (auto i : dimensionRange())
        dWdLambda[i] = (lambdaBar[i] * dWdLambdaBar[i] - (1.0 / 3.0) * sumLambdaBar) / lambda[i];

    } else {
      for (auto j : parameterRange())
        for (auto k : dimensionRange())
          dWdLambda[k] += (mu[j] * (pow(lambda[k], alpha[j]) - 1)) / lambda[k];
    }
    return dWdLambda;
  }

  /**
   * \brief Computes the second derivatives of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  SecondDerivative secondDerivativeImpl(const PrincipalStretches& lambda) const {
    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;
    auto dS     = SecondDerivative::Zero().eval();

    if constexpr (usesDeviatoricStretches) {
      const auto lambdaBar = Impl::deviatoricStretches(lambda);
      const auto dWdLambda = firstDerivativeImpl(lambda);

      for (auto a : dimensionRange()) {
        for (auto b : dimensionRange()) {
          if (a == b) {
            for (auto p : parameterRange()) {
              ScalarType sumC{0.0};
              for (auto c : dimensionRange())
                sumC += pow(lambdaBar[c], alpha[p]);
              dS(a, b) += mu[p] * alpha[p] * ((1.0 / 3.0) * pow(lambdaBar[a], alpha[p]) + (1.0 / 9.0) * sumC);
            }
          } else {
            for (auto p : parameterRange()) {
              ScalarType sumC{0.0};
              for (auto c : dimensionRange())
                sumC += pow(lambdaBar[c], alpha[p]);
              dS(a, b) +=
                  mu[p] * alpha[p] *
                  (-(1.0 / 3.0) * (pow(lambdaBar[a], alpha[p]) + pow(lambdaBar[b], alpha[p])) + (1.0 / 9.0) * sumC);
            }
          }

          dS(a, b) *= 1.0 / (lambda[a] * lambda[b]);

          if (a == b)
            dS(a, b) -= (2.0 / lambda[a]) * dWdLambda[a];
        }
      }
    } else {
      for (auto j : parameterRange())
        for (auto k : dimensionRange())
          dS(k, k) += (-2 * (mu[j] * (pow(lambda[k], alpha[j]) - 1)) + (mu[j] * pow(lambda[k], alpha[j]) * alpha[j])) /
                      pow(lambda[k], 2);
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
    return OgdenT<STO, nOgdenParameters, stretchTag>(materialParameters_, ogdenParameters_);
  }

private:
  MaterialParameters materialParameters_;
  OgdenParameters ogdenParameters_;

  inline auto parameterRange() const { return Dune::Hybrid::integralRange(nOgdenParameters); }
  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(dim); }
};

/**
 * \brief Alias for OgdenT with double as the default scalar type.
 */
template <int n, PrincipalStretchTag tag>
using Ogden = OgdenT<double, n, tag>;

} // namespace Ikarus::Materials
