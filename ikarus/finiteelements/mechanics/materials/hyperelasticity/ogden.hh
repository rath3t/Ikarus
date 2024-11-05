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

  ScalarType storedEnergyImpl(const PrincipalStretches& lambda) const {
    auto lambdaBar = transformStretches(lambda);

    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;

    ScalarType energy{};
    for (auto i : parameterRange())
      energy += mu[i] / alpha[i] *
                (pow(lambdaBar[0], alpha[i]) + pow(lambdaBar[1], alpha[i]) + pow(lambdaBar[2], alpha[i]) - 3);

    return energy;
  }

  FirstDerivative firstDerivativeImpl(const PrincipalStretches& lambda) const {
    auto lambdaBar = transformStretches(lambda);

    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;

    auto dWdLambdaBar = FirstDerivative::Zero().eval();
    for (auto j : parameterRange())
      for (auto k : dimensionRange())
        dWdLambdaBar[k] += mu[j] * (pow(lambdaBar[k], alpha[j] - 1));

    ScalarType sumLambdaBar{0.0};
    for (auto b : dimensionRange())
      sumLambdaBar += lambdaBar[b] * dWdLambdaBar[b];

    auto dWdLambda = FirstDerivative::Zero().eval();
    for (auto i : dimensionRange())
      dWdLambda[i] = (lambdaBar[i] * dWdLambdaBar[i] - (1.0 / 3.0) * sumLambdaBar) / lambda[i];

    return dWdLambda;
  }

  SecondDerivative secondDerivativeImpl(const PrincipalStretches& lambda) const {
    auto lambdaBar = transformStretches(lambda);
    auto S         = firstDerivativeImpl(lambda);

    auto& mu    = materialParameters_;
    auto& alpha = ogdenParameters_;

    // for (auto j : parameterRange())
    //   for (auto k : dimensionRange())
    //     dS(k, k) += mu[j] * alpha[j]  /* * (alpha[j] - 1)*/ * pow(lambdas[k], alpha[j] - 2);

    auto dS = SecondDerivative::Zero().eval();

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

    return dS;
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return OgdenT<ScalarTypeOther> The rebound Ogden material.
   */
  template <typename STO>
  auto rebind() const {
    return OgdenT<STO, nOgdenParameters>(materialParameters_, ogdenParameters_);
  }

private:
  MaterialParameters materialParameters_;
  OgdenParameters ogdenParameters_;

  inline auto parameterRange() const { return Dune::Hybrid::integralRange(nOgdenParameters); }
  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(worldDimension); }

  PrincipalStretches transformStretches(const PrincipalStretches& lambda) const {
    ScalarType J    = std::accumulate(lambda.begin(), lambda.end(), ScalarType{1.0}, std::multiplies());
    ScalarType Jmod = pow(J, -1.0 / 3.0);

    auto lambdaBar = PrincipalStretches::Zero().eval();
    for (auto i : dimensionRange())
      lambdaBar[i] = Jmod * lambda[i];

    return lambdaBar;
  }
};

/**
 * \brief Alias for OgdenT with double as the default scalar type.
 */
template <int n>
using Ogden = OgdenT<double, n>;

} // namespace Ikarus
