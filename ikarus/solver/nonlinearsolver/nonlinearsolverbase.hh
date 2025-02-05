// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverbase.hh
 * \brief Base for all nonlinear solvers
 */

#pragma once
#include <type_traits>

#include <ikarus/solver/nonlinearsolver/nonlinearsolverstate.hh>
#include <ikarus/utils/broadcaster/broadcaster.hh>
#include <ikarus/utils/broadcaster/broadcastermessages.hh>

namespace Ikarus {
// namespace Impl {

#define COMMON_SIGNATURES                                                                                   \
  void(NonLinearSolverMessages), void(NonLinearSolverMessages, double), void(NonLinearSolverMessages, int), \
      void(NonLinearSolverMessages, typename NonlinearSolverStateFactory<NLO>::type&)

template <typename NLO>
struct NonlinearSolverBase : public Broadcasters<COMMON_SIGNATURES>
{
};

} // namespace Ikarus
