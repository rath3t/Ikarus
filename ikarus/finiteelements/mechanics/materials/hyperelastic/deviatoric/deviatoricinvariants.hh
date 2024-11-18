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
    ScalarType W1                = invariants[0] * pow(invariants[2], -1.0 / 3.0);
    ScalarType W2                = invariants[1] * pow(invariants[2], -2.0 / 3.0);
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

} // namespace Ikarus::Materials