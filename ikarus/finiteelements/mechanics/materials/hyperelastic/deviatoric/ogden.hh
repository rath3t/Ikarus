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
 * \ingroup materials
 *
 * \details If total principal stretches are used, the energy is computed via
 * \f[ \hat{\Psi}(\lambda_1, \lambda_2, \lambda_3) = \sum_{p=1}^n{
 \frac{\mu_n}{\alpha_n} ({\lambda_1}^{\alpha_n} + {\lambda_2}^{\alpha_n} +
 {\lambda_3}^{\alpha_n} - 3) - \mu_n \ln J}, \f]
 * with \f$ J = \lambda_1 \lambda_2 \lambda_3 \f$.
 *
 * If deviatoric principal stretches (\f$ \bar{\lambda_i} = \lambda_i J^{-1/3} \f$)
 * are used, the energy is computed via
 * \f[ \hat{\Psi}(\lambda_1, \lambda_2, \lambda_3) = \sum_{p=1}^n{
 \frac{\mu_n}{\alpha_n}({\bar{\lambda_1}}^{\alpha_n} + {\bar{\lambda_2}}^{\alpha_n} +
 {\bar{\lambda_3}}^{\alpha_n} - 3) }. \f]
 *
 * \remark See \cite bergstromMechanicsSolidPolymers2015 for details on this material.
 * \tparam ST The underlying scalar type.
 * \tparam n Number of ogden parameters
 * \tparam tag Type of principal stretch quantity, either total stretches or deviatoric stretches
 */
template <typename ST, int n, PrincipalStretchTag tag>
struct OgdenT
{
  using ScalarType         = ST;
  using PrincipalStretches = Eigen::Vector<ScalarType, 3>;

  static constexpr PrincipalStretchTag stretchTag = tag;
  static constexpr int numMatParameters           = n;
  static constexpr int dim                        = 3;
  static constexpr bool usesDeviatoricStretches   = stretchTag == PrincipalStretchTag::deviatoric;

  using FirstDerivative  = Eigen::Vector<ScalarType, dim>;
  using SecondDerivative = Eigen::Matrix<ScalarType, dim, dim>;

  using MaterialParameters = std::array<double, numMatParameters>;
  using MaterialExponents  = std::array<double, numMatParameters>;

  [[nodiscard]] constexpr static std::string name() noexcept {
    return "Ogden (n = " + std::to_string(numMatParameters) + ", stretch type = " + toString(tag) + ")";
  }

  /**
   * \brief Constructor for OgdenT
   *
   * \param mpt material parameters (array of mu values)
   * \param mex material exponents (array of alpha values)
   */
  explicit OgdenT(const MaterialParameters& mpt, const MaterialExponents& mex)
      : materialParameters_{mpt},
        materialExponents_{mex} {}

  /**
   * \brief Returns the material parameters (mu values) stored in the material
   */
  const MaterialParameters& materialParametersImpl() const { return materialParameters_; }

  /**
   * \brief Returns the material exponents (alpha values) stored in the material
   */
  const MaterialExponents& materialExponents() const { return materialExponents_; }

  /**
   * \brief Computes the stored energy in the Ogden material model.
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  ScalarType storedEnergyImpl(const PrincipalStretches& lambda) const {
    auto& mu    = materialParameters_;
    auto& alpha = materialExponents_;

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
    auto& alpha    = materialExponents_;
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
    auto& alpha = materialExponents_;
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
    return OgdenT<STO, numMatParameters, stretchTag>(materialParameters_, materialExponents_);
  }

private:
  MaterialParameters materialParameters_;
  MaterialExponents materialExponents_;

  inline auto parameterRange() const { return Dune::Hybrid::integralRange(numMatParameters); }
  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(dim); }
};

/**
 * \brief Alias for OgdenT with double as the default scalar type.
 */
template <int n, PrincipalStretchTag tag>
using Ogden = OgdenT<double, n, tag>;

} // namespace Ikarus::Materials
