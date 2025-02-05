// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverstate.hh
 * \brief State for all nonlinear solvers
 */

#pragma once

#include <Eigen/Core>

namespace Ikarus {

template <typename P1, typename P2, typename SolType = Eigen::VectorXd>
struct NonlinearSolverState
{
  const P1 firstParameter;  // Solution
  const P2 secondParameter; //
  const SolType solution;

  double rNorm{};
  double dNorm{};
  int iteration{};
};

template <typename NLO>
struct NonlinearSolverStateFactory
{
private:
  static constexpr bool derivativeTypeEqualValueType =
      std::is_same_v<typename NLO::ValueType, typename NLO::DerivativeType>;
  // using SolutionType = std::conditional_t<derivativeTypeEqualValueType, const typename NLO::ValueType&,
  //                                         const typename NLO::DerivativeType&>;
  using SolutionType = const typename NLO::ValueType&;

public:
  using type = NonlinearSolverState<const typename NLO::template ParameterValue<0>&,
                                    typename NLO::template ParameterValue<1>, SolutionType>;
};

} // namespace Ikarus
