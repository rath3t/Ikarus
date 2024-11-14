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
    const auto& pex = exponents(pex_);
    const auto& qex = exponents(qex_);

    for (auto i : parameterRange())
      energy += matParameters_[i] * pow((invariants[0] * pow(invariants[2], -1.0 / 3.0) - 3.0), pex[i]) *
                pow((invariants[1] * pow(invariants[2], -2.0 / 3.0) - 3.0), qex[i]);

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

    const auto& mu  = matParameters_;
    const auto& pex = exponents(pex_);
    const auto& qex = exponents(qex_);

    ScalarType I1        = invariants[0];
    ScalarType I2        = invariants[1];
    ScalarType I3        = invariants[2];
    ScalarType I3Pow1by3 = pow(I3, 1.0 / 3.0);
    ScalarType I3Pow2by3 = pow(I3, 2.0 / 3.0);
    ScalarType I3Pow4by3 = pow(I3, 4.0 / 3.0);
    ScalarType I4        = -I1 + 3.0 * I3Pow1by3;
    ScalarType I5        = -I2 + 3.0 * I3Pow2by3;

    for (auto j : parameterRange())
      for (auto k : dimensionRange()) {
        auto factor1 = -2.0 * mu[j] * pow(-I4 / I3Pow1by3, pex[j]) * pow(-I5 / I3Pow2by3, qex[j]);
        auto factor2 = 1.0 / (pow(lambda[k], 3.0) * I4 * I5);
        auto factor3 = (-1.0 / 3.0) * pow(lambda[k], 2.0) * I5 * (I1 - 3.0 * pow(lambda[k], 2.0));
        auto factor4 = I1 * I3 - 3.0 * I3Pow4by3 + I2 * I3Pow1by3 * pow(lambda[k], 2.0) -
                       (1.0 / 3.0) * I1 * I2 * pow(lambda[k], 2.0);
        dWdLambda[k] += factor1 * factor2 * (factor3 * pex[j] + factor4 * qex[j]);
      }

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
