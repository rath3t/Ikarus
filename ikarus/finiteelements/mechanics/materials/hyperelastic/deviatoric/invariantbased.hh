// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file ogden.hh
 * \brief Implementation of the regularized InvariantBased material model.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

namespace Impl {
  template <typename PS>
  struct DeviatoricInvariants
  {
    using PrincipalStretches = PS;
    using ScalarType         = PrincipalStretches::value_type;
    static constexpr int dim = PrincipalStretches::RowsAtCompileTime;
    using Invariants         = PrincipalStretches;

    using FirstDerivative  = Eigen::Vector<ScalarType, dim>;
    using SecondDerivative = Eigen::Matrix<ScalarType, dim, dim>;

    explicit DeviatoricInvariants(const PrincipalStretches& lambda)
        : lambda_{lambda} {}

    auto value() const {
      const Invariants& invariants = Impl::invariants(lambda_);
      ScalarType W1                = invariants[0] * pow(invariants[2], -1.0 / 3.0) - 3.0;
      ScalarType W2                = invariants[1] * pow(invariants[2], -2.0 / 3.0) - 3.0;
      return std::make_pair(W1, W2);
    }

    auto firstDerivative() const {
      const Invariants& invariants            = Impl::invariants(lambda_);
      auto dW1dLambda                         = FirstDerivative::Zero().eval();
      auto dW2dLambda                         = FirstDerivative::Zero().eval();
      auto [I1, I2, I3, I3Pow1by3, I3Pow2by3] = computeInvariants(invariants);

      for (auto i : dimensionRange()) {
        dW1dLambda[i] = 2.0 * (3.0 * pow(lambda_[i], 2.0) - I1) / (3.0 * lambda_[i] * I3Pow1by3);
        dW2dLambda[i] = -2.0 * (3.0 * (I3 / pow(lambda_[i], 2.0)) - I2) / (3.0 * lambda_[i] * I3Pow2by3);
      }

      return std::make_pair(dW1dLambda, dW2dLambda);
    }

    auto secondDerivative() const {
      const Invariants& invariants            = Impl::invariants(lambda_);
      auto ddW1dLambda                        = SecondDerivative::Zero().eval();
      auto ddW2dLambda                        = SecondDerivative::Zero().eval();
      auto [I1, I2, I3, I3Pow1by3, I3Pow2by3] = computeInvariants(invariants);

      for (auto i : dimensionRange())
        for (auto j : dimensionRange()) {
          if (i == j) {
            ddW1dLambda(i, j) = (2.0 / 9.0) * (5 * I1 - 3 * pow(lambda_[i], 2.0)) / (pow(lambda_[i], 2.0) * I3Pow1by3);
            ddW2dLambda(i, j) =
                (2.0 / 9.0) * ((15.0 * I3 / pow(lambda_[i], 2.0)) - I2) / (pow(lambda_[i], 2.0) * I3Pow2by3);
          } else {
            ddW1dLambda(i, j) = (4.0 / 9.0) * (I1 - 3 * (pow(lambda_[i], 2.0) + pow(lambda_[j], 2.0))) /
                                (lambda_[i] * lambda_[j] * I3Pow1by3);
            ddW2dLambda(i, j) = (-4.0 / 9.0) * (2.0 * I2 - 3 * pow(lambda_[i], 2.0) * pow(lambda_[j], 2.0)) /
                                (lambda_[i] * lambda_[j] * I3Pow2by3);
          }
        }
      return std::make_pair(ddW1dLambda, ddW2dLambda);
    }

  private:
    PrincipalStretches lambda_;

    inline auto dimensionRange() const { return Dune::Hybrid::integralRange(dim); }

    auto computeInvariants(const Invariants& invariants) const {
      ScalarType I1        = invariants[0];
      ScalarType I2        = invariants[1];
      ScalarType I3        = invariants[2];
      ScalarType I3Pow1by3 = pow(I3, 1.0 / 3.0);
      ScalarType I3Pow2by3 = pow(I3, 2.0 / 3.0);

      return std::make_tuple(I1, I2, I3, I3Pow1by3, I3Pow2by3);
    }
  };

} // namespace Impl

/**
 * \brief Implementation of the InvariantBased material model.
 *
 * \tparam ST The scalar type for the strains and stresses,....
 * \tparam n number of ogden parameters
 * \tparam tag type of principal stretch quantity, either total stretches or deviatoric stretches
 * \ingroup materials
 */
template <typename ST, int n>
struct InvariantBasedT
{
  static constexpr int dim              = 3;
  using ScalarType                      = ST;
  using PrincipalStretches              = Eigen::Vector<ScalarType, dim>;
  using Invariants                      = PrincipalStretches;
  static constexpr int numMatParameters = n;

  using Exponents          = std::array<std::size_t, numMatParameters>;
  using MaterialParameters = std::array<double, numMatParameters>;

  using FirstDerivative  = Eigen::Vector<ScalarType, dim>;
  using SecondDerivative = Eigen::Matrix<ScalarType, dim, dim>;

  [[nodiscard]] constexpr static std::string name() noexcept {
    return "InvariantBased (n = " + std::to_string(numMatParameters);
  }

  const MaterialParameters& materialParametersImpl() const { return matParameters_; }

  const Exponents& pExponents() const { return pex_; }

  const Exponents& qExponents() const { return qex_; }

  /**
   * \brief Constructor for InvariantBasedT
   *
   * \param mpt material parameters (array of mu values)
   * \param opt ogden parameters (array of alpha values)
   */
  explicit InvariantBasedT(const Exponents& pex, const Exponents& qex, const MaterialParameters& matParameters)
      : pex_{pex},
        qex_{qex},
        matParameters_{matParameters} {}

  /**
   * \brief Computes the stored energy in the InvariantBased material model.
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  ScalarType storedEnergyImpl(const PrincipalStretches& lambda) const {
    const Invariants& invariants = Impl::invariants(lambda);
    ScalarType energy{0.0};
    const auto& pex           = exponents(pex_);
    const auto& qex           = exponents(qex_);
    const auto& devInvariants = Impl::DeviatoricInvariants<PrincipalStretches>(lambda);
    const auto& [W1, W2]      = devInvariants.value();

    for (auto i : parameterRange())
      energy += matParameters_[i] * pow(W1, pex[i]) * pow(W2, qex[i]);

    return energy;
  }

  /**
   * \brief Computes the first derivative of the stored energy function w.r.t. the total principal stretches
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  FirstDerivative firstDerivativeImpl(const PrincipalStretches& lambda) const {
    const Invariants& invariants = Impl::invariants(lambda);
    auto dWdLambda               = FirstDerivative::Zero().eval();

    const auto& mu                       = matParameters_;
    const auto& pex                      = exponents(pex_);
    const auto& qex                      = exponents(qex_);
    const auto& devInvariants            = Impl::DeviatoricInvariants<PrincipalStretches>(lambda);
    const auto& [W1, W2]                 = devInvariants.value();
    const auto& [dW1dLambda, dW2dLambda] = devInvariants.firstDerivative();

    for (auto j : parameterRange())
      for (auto k : dimensionRange())
        dWdLambda[k] +=
            mu[j] * pow(W1, pex[j]) * pow(W2, qex[j]) * ((pex[j] * dW1dLambda[k]) / W1 + (qex[j] * dW2dLambda[k]) / W2);

    return dWdLambda;
  }

  /**
   * \brief Computes the second derivatives of the stored energy function w.r.t. the total principal stretches
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  SecondDerivative secondDerivativeImpl(const PrincipalStretches& lambda) const {
    const Invariants& invariants = Impl::invariants(lambda);
    auto dS                      = SecondDerivative::Zero().eval();

    const auto& mu                         = matParameters_;
    const auto& pex                        = exponents(pex_);
    const auto& qex                        = exponents(qex_);
    const auto& devInvariants              = Impl::DeviatoricInvariants<PrincipalStretches>(lambda);
    const auto& [W1, W2]                   = devInvariants.value();
    const auto& [dW1dLambda, dW2dLambda]   = devInvariants.firstDerivative();
    const auto& [ddW1dLambda, ddW2dLambda] = devInvariants.secondDerivative();

    for (auto p : parameterRange())
      for (auto i : dimensionRange())
        for (auto j : dimensionRange()) {
          auto factor1 = (pex[p] / W1) * (((pex[p] - 1.0) / W1) * dW1dLambda[i] * dW1dLambda[j] + ddW1dLambda(i, j));
          auto factor2 = (qex[p] / W2) * (((qex[p] - 1.0) / W2) * dW2dLambda[i] * dW2dLambda[j] + ddW2dLambda(i, j));
          auto factor3 =
              ((pex[p] * qex[p]) / (W1 * W2)) * (dW1dLambda[i] * dW2dLambda[j] + dW1dLambda[j] * dW2dLambda[i]);
          auto factor4 = mu[p] * pow(W1, pex[p]) * pow(W2, qex[p]);
          dS(i, j) += (1.0 / lambda[i]) * (factor4 * (factor1 + factor2 + factor3));
          if (i == j) {
            auto factor5 = (pex[p] / W1 * dW1dLambda[i]) + (qex[p] / W2 * dW2dLambda[i]);
            dS(i, j) -= (1.0 / pow(lambda[i], 2.0)) * factor4 * factor5;
          }
        }
    return dS;
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return InvariantBasedT<ScalarTypeOther> The rebound InvariantBased material.
   */
  template <typename STO>
  auto rebind() const {
    return InvariantBasedT<STO, numMatParameters>(pex_, qex_, matParameters_);
  }

private:
  Exponents pex_, qex_;
  MaterialParameters matParameters_;

  inline auto parameterRange() const { return Dune::Hybrid::integralRange(numMatParameters); }
  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(dim); }

  /** \brief A function to transform the underlying type of the exponents to ScalarType.
   *
   * \details If Concepts::AutodiffScalar<ScalarType> is true, then the pow function returns zero instead of one, if
   * a certain base is raised to an exponent of zero. Hence this transformation is necessary.
   * \param exp_ An array of exponents of type std::size_t.
   * \return An array of exponents of type \tparam ScalarType if Concepts::AutodiffScalar<ScalarType> is true, else the
   * given \param exp_ is returned..
   */
  decltype(auto) exponents(const Exponents& exp_) const {
    if constexpr (Concepts::AutodiffScalar<ScalarType>) {
      std::array<ScalarType, numMatParameters> transformedExp_{};
      std::ranges::transform(exp_, transformedExp_.begin(), [](std::size_t x) { return static_cast<ScalarType>(x); });
      return transformedExp_;
    } else
      return exp_;
  }
};

/**
 * \brief Alias for InvariantBasedT with double as the default scalar type.
 */
template <int n>
using InvariantBased = InvariantBasedT<double, n>;

} // namespace Ikarus::Materials
