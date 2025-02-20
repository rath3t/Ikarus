// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverbase.hh
 * \brief Base for all nonlinear solvers
 */

#pragma once
#include <ikarus/solver/nonlinearsolver/nonlinearsolverstate.hh>
#include <ikarus/utils/broadcaster/broadcaster.hh>
#include <ikarus/utils/broadcaster/broadcastermessages.hh>

namespace Ikarus {

/**
 * \brief Base for all nonlinear solvers. Defines the message interface that can be broadcasted to listeners.
 *
 * \tparam NLO The nonlinear operator
 * \tparam Args Additional message signatures can be broadcasted
 */
template <typename NLO, typename... Args>
struct NonlinearSolverBase : public Broadcasters<void(NonLinearSolverMessages), void(NonLinearSolverMessages, double),
                                                 void(NonLinearSolverMessages, int),
                                                 void(NonLinearSolverMessages, NonlinearSolverStateType<NLO>&), Args...>
{
  using State = NonlinearSolverStateType<NLO>;
};

} // namespace Ikarus
