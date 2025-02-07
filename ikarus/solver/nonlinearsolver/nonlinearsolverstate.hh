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
template <typename TypeListOne, typename TypeListTwo>
class NonLinearOperator;

enum class FESolutions;
enum class FEParameter;
template <FESolutions sol, FEParameter para, typename SV, typename PM>
class FERequirements;

template <typename P1, typename SolType = P1>
struct NonlinearSolverState
{
  P1 correction;
  SolType solution;

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
    // For NLOs which have more than 2 functions, we use the derivative as P1 (e.g. TR)
    static constexpr bool useDerivativeType = std::tuple_size_v<typename NLO::FunctionReturnValues> > 2;
    using P1Type =
        std::conditional_t<useDerivativeType, const typename NLO::DerivativeType&, const typename NLO::ValueType&>;

  public:
    using type = NonlinearSolverState<P1Type, const typename NLO::template ParameterValue<0>&>;
  };
} // namespace Impl

template <typename T>
using NonlinearSolverStateType = Impl::NonlinearSolverStateFactory<T>::type;

} // namespace Ikarus
