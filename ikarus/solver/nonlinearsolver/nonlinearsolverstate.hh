// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverstate.hh
 * \brief State for all nonlinear solvers
 */

#pragma once

#include <ikarus/utils/traits.hh>

namespace Ikarus {

// Forward declarations
#ifndef DOXYGEN
template <typename TypeListOne, typename TypeListTwo>
class NonLinearOperator;
#endif

/**
 * \brief State for nonlinear solvers
 *
 * \tparam Correction the type of the correction vector
 * \tparam SolutionType the type of the solution vector
 */
template <typename Correction, typename SolutionType = Correction>
struct NonlinearSolverState
{
  Correction correction;
  SolutionType solution;

  double rNorm{};
  double dNorm{};
  int iteration{};
};

namespace Impl {
  template <typename T>
  struct NonlinearSolverStateFactory;

  template <typename NLO>
  requires traits::isSpecialization<NonLinearOperator, NLO>::value
  struct NonlinearSolverStateFactory<NLO>
  {
  private:
    // For NLOs which have more than 2 functions, we use the derivative as Correction (e.g. TR)
    static constexpr bool useDerivativeType = std::tuple_size_v<typename NLO::FunctionReturnValues> > 2;
    using CorrectionType =
        std::conditional_t<useDerivativeType, const typename NLO::DerivativeType&, const typename NLO::ValueType&>;
    using SolutionType = const typename NLO::template ParameterValue<0>&;

  public:
    using type = NonlinearSolverState<CorrectionType, SolutionType>;
  };
} // namespace Impl

/**
 * \brief Helper to deduce the correct types for NonlinearSolverState
 *
 * \tparam NLO The nonlinear operator
 */
template <typename NLO>
using NonlinearSolverStateType = Impl::NonlinearSolverStateFactory<NLO>::type;

} // namespace Ikarus
