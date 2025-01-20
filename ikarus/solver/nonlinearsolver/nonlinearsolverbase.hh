// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverbase.hh
 * \brief Base for all nonlinear solvers
 */

#pragma once
#include <type_traits>

#include <ikarus/utils/broadcaster/broadcaster.hh>
#include <ikarus/utils/broadcaster/broadcastermessages.hh>

namespace Ikarus {
namespace Impl {

#define COMMON_SIGNATURES                                                                                   \
  void(NonLinearSolverMessages, typename NLO::template ParameterValue<0>&, const typename NLO::ValueType&), \
      void(NonLinearSolverMessages), void(NonLinearSolverMessages, double), void(NonLinearSolverMessages, int)

  // Helper struct to conditionally add the second signature
  template <typename NLO, bool IsSame>
  struct SignatureHelper;

  // Specialization when ValueType and DerivativeType are the same
  template <typename NLO>
  struct SignatureHelper<NLO, true>
  {
    using Type = Broadcasters<COMMON_SIGNATURES>;
  };

  // Specialization when ValueType and DerivativeType are different
  template <typename NLO>
  struct SignatureHelper<NLO, false>
  {
    using Type = Broadcasters<void(NonLinearSolverMessages, typename NLO::template ParameterValue<0>&,
                                   const typename NLO::DerivativeType&),
                              COMMON_SIGNATURES>;
  };
} // namespace Impl
template <typename NLO>
struct NonlinearSolverBase
    : public Impl::SignatureHelper<NLO,
                                   std::is_same<typename NLO::ValueType, typename NLO::DerivativeType>::value>::Type
{
};

} // namespace Ikarus
