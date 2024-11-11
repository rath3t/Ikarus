// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file regularizedhelpers.hh
 * \brief Implementation of helper functions for regularized material models.
 * \ingroup  materials
 */

#pragma once

#include <functional>
#include <numeric>

#include <dune/common/hybridutilities.hh>

namespace Ikarus::Impl {

template <typename PrincipalStretches>
inline PrincipalStretches regularizeStretches(const PrincipalStretches& lambda) {
  using ScalarType = std::remove_cvref_t<decltype(lambda[0])>;

  ScalarType J    = std::accumulate(lambda.begin(), lambda.end(), ScalarType{1.0}, std::multiplies());
  ScalarType Jmod = pow(J, -1.0 / 3.0);

  auto lambdaBar = PrincipalStretches::Zero().eval();
  for (auto i : Dune::Hybrid::integralRange(3))
    lambdaBar[i] = Jmod * lambda[i];

  return lambdaBar;
}

} // namespace Ikarus::Impl