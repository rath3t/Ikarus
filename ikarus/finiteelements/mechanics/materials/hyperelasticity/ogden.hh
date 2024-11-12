// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file ogden.hh
 * \brief Implementation of the regularized Ogden material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

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
    return "Ogden (n=" + std::to_string(nOgdenParameters) + ")";
  }

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

  SecondDerivative secondDerivativeImpl(const PrincipalStretches& lambda) const {
    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;
    auto dS     = SecondDerivative::Zero().eval();

    if constexpr (usesDeviatoricStretches) {
      auto lambdaBar = Impl::deviatoricStretches(lambda);
      auto S         = firstDerivativeImpl(lambda);

      for (auto a : dimensionRange())
        for (auto b : dimensionRange()) {
          if (a == b) {
            for (auto p : parameterRange()) {
              ScalarType sumC{0.0};
              for (auto c : dimensionRange())
                sumC += pow(lambdaBar[c], alpha[p]);

              dS(a, b) = mu[p] * alpha[p] * ((1.0 / 3.0) * pow(lambdaBar[a], alpha[p]) + (1.0 / 9.0) * sumC);
            }
          } else {
            for (auto p : parameterRange()) {
              ScalarType sumC{0.0};
              for (auto c : dimensionRange())
                sumC += pow(lambdaBar[c], alpha[p]);

              dS(a, b) = mu[p] * alpha[p] *
                         (-(1.0 / 3.0) * pow(lambdaBar[a], alpha[p]) - (1.0 / 3.0) * pow(lambdaBar[b], alpha[p]) +
                          (1.0 / 9.0) * sumC);
            }
          }

          auto factor = pow(lambda[a], -2) / lambda[b];
          dS(a, b) *= factor;

          if (a == b) {
            ScalarType sumSiso{0.0};
            for (auto i : dimensionRange())
              sumSiso += 2 / pow(lambda[i], 2) * (S[i] / lambda[i]);

            dS(a, b) -= sumSiso;
          }
        }
    } else {
      for (auto j : parameterRange())
        for (auto k : dimensionRange())
          dS(k, k) += -2 * (mu[j] * (pow(lambda[k], alpha[j]) - 1)) / pow(lambda[k], 3) +
                      (mu[j] * pow(lambda[k], alpha[j]) * alpha[j] / lambda[k]) / pow(lambda[k], 2);
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
