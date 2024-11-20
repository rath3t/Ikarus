// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file invariantbased.hh
 * \brief Implementation of the InvariantBased material model.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/deviatoricinvariants.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Implementation of the InvariantBased material model.
 * \ingroup materials
 *
 * \details This is a general model based on deviatoric invariants. It can be used to derive other specific material
 * models, for instance Mooney-Rivlin and Yeoh models.
 * The energy is computed as
 * \f[ \hat{\Psi}(\lambda_1, \lambda_2, \lambda_3) = \sum_{p,q=0}^n{
 C_{pq} (W_1 - 3)^p (W_2 - 3)^q}. \f]
 *
 * \remark See \cite bergstromMechanicsSolidPolymers2015 for details on this material. For information on the deviatoric
 * invariant \f$ W_1 \f$, see \ref DeviatoricInvariants.
 *
 * \tparam ST The underlying scalar type.
 * \tparam n Number of material parameters
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

  static constexpr auto stretchTag = PrincipalStretchTag::deviatoric;

  [[nodiscard]] constexpr static std::string name() noexcept {
    return "InvariantBased (n = " + std::to_string(numMatParameters);
  }

  const MaterialParameters& materialParametersImpl() const { return matParameters_; }

  const Exponents& pExponents() const { return pex_; }

  const Exponents& qExponents() const { return qex_; }

  /**
   * \brief Constructor for InvariantBasedT
   *
   * \param pex Array of exponents related to the first invariant
   * \param qex Array of exponents related to the second invariant
   * \param matParameters Array of material parameters
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
    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches>(lambda);
    auto [W1, W2]             = devInvariants.value();
    W1 -= 3.0;
    W2 -= 3.0;

    for (auto i : parameterRange())
      energy += matParameters_[i] * pow(W1, pex[i]) * pow(W2, qex[i]);

    return energy;
  }

  /**
   * \brief Computes the first derivative of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  FirstDerivative firstDerivativeImpl(const PrincipalStretches& lambda) const {
    const Invariants& invariants = Impl::invariants(lambda);
    auto dWdLambda               = FirstDerivative::Zero().eval();

    const auto& mu            = matParameters_;
    const auto& pex           = exponents(pex_);
    const auto& qex           = exponents(qex_);
    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches>(lambda);
    auto [W1, W2]             = devInvariants.value();
    W1 -= 3.0;
    W2 -= 3.0;
    const auto& [dW1dLambda, dW2dLambda] = devInvariants.firstDerivative();

    for (auto p : parameterRange())
      for (auto k : dimensionRange()) {
        auto W1pm1p = safeMultiply(pow(W1, pex[p] - 1.0), pex[p]);
        auto W2qm1q = safeMultiply(pow(W2, qex[p] - 1.0), qex[p]);
        dWdLambda[k] +=
            mu[p] * ((W1pm1p * pow(W2, qex[p]) * dW1dLambda[k]) + (pow(W1, pex[p]) * W2qm1q * dW2dLambda[k]));
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
    const Invariants& invariants = Impl::invariants(lambda);
    auto dS                      = SecondDerivative::Zero().eval();

    const auto& mu            = matParameters_;
    const auto& pex           = exponents(pex_);
    const auto& qex           = exponents(qex_);
    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches>(lambda);
    auto [W1, W2]             = devInvariants.value();
    W1 -= 3.0;
    W2 -= 3.0;
    const auto& [dW1dLambda, dW2dLambda]   = devInvariants.firstDerivative();
    const auto& [ddW1dLambda, ddW2dLambda] = devInvariants.secondDerivative();

    for (auto p : parameterRange())
      for (auto i : dimensionRange())
        for (auto j : dimensionRange()) {
          auto W1pm1p       = safeMultiply(pow(W1, pex[p] - 1.0), pex[p]);
          auto W2qm1q       = safeMultiply(pow(W2, qex[p] - 1.0), qex[p]);
          auto W1pm2pp      = safeMultiply(pow(W1, pex[p] - 2.0), pex[p] * (pex[p] - 1.0));
          auto W2qm2qq      = safeMultiply(pow(W2, qex[p] - 2.0), qex[p] * (qex[p] - 1.0));
          auto dW1W2dlambda = dW1dLambda[i] * dW2dLambda[j] + dW1dLambda[j] * dW2dLambda[i];
          auto factor1      = (W2qm1q * dW1W2dlambda + pow(W2, qex[p]) * ddW1dLambda(i, j)) * W1pm1p * mu[p];
          auto factor2      = W2qm1q * ddW2dLambda(i, j) * pow(W1, pex[p]) * mu[p];
          auto factor3      = W1pm2pp * pow(W2, qex[p]) * dW1dLambda[i] * dW1dLambda[j] * mu[p];
          auto factor4      = W2qm2qq * pow(W1, pex[p]) * dW2dLambda[i] * dW2dLambda[j] * mu[p];
          dS(i, j) += factor1 + factor2 + factor3 + factor4;
          if (i == j) {
            auto factor5 =
                mu[p] * (W1pm1p * pow(W2, qex[p]) * dW1dLambda[i] + pow(W1, pex[p]) * W2qm1q * dW2dLambda[i]);
            dS(i, j) -= (1.0 / lambda[i]) * factor5;
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

  /**
   * \brief A function to safely multiply infinity times zero and return zero instead of nan.
   *
   * \details In the undeformed configuration, all principal stretches are equal to one. In this scenario, the
   * deviatoric invariants W1 and W2 are zero. For any positive integer n, 0^{-n} is inf. Multiplying this with 0.0
   * again leads to nan. In order to circumvent this, in such a scenario, zero is returned instead of nan.
   *
   * \remark If inf is multiplied with any other number apart from zero, nan is returned.
   *
   * \param x First number to be multiplied.
   * \param y Second number to be multiplied.
   * \return ScalarType Either 0.0 or x * y.
   */
  ScalarType safeMultiply(ScalarType x, ScalarType y) const {
    auto checkInfinityTimeZero = [](double val1, double val2) -> bool {
      constexpr double tol = 1e-14;
      return ((std::isinf(val1) && Dune::FloatCmp::eq(val2, 0.0, tol)) ||
              (std::isinf(val2) && Dune::FloatCmp::eq(val1, 0.0, tol)));
    };

    // static_cast needed if Concepts::AutodiffScalar<ScalarType>
    return checkInfinityTimeZero(static_cast<double>(x), static_cast<double>(y)) ? ScalarType{0.0} : x * y;
  }
};

/**
 * \brief Alias for InvariantBasedT with double as the default scalar type.
 */
template <int n>
using InvariantBased = InvariantBasedT<double, n>;

} // namespace Ikarus::Materials
