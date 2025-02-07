// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controlroutinebase.hh
 * \brief Base for all control routines
 */

#pragma once
#include <ikarus/controlroutines/controlroutinestate.hh>
#include <ikarus/utils/broadcaster/broadcaster.hh>
#include <ikarus/utils/broadcaster/broadcastermessages.hh>

namespace Ikarus {

/**
 * \brief Base for all control routines. Defines the message interface that can be broadcasted to listeners.
 *
 * \tparam NLO The nonlinear operator
 * \tparam Args Additional message signatures can be broadcasted
 */
template <typename NLO, typename... Args>
struct ControlRoutineBase
    : public Broadcasters<void(ControlMessages), void(ControlMessages, const std::string&),
                          void(ControlMessages, int, const std::string&), void(ControlMessages, int, double),
                          void(ControlMessages, ControlRoutineStateType<NLO>&), Args...>

{
  using State = ControlRoutineStateType<NLO>;
};

} // namespace Ikarus
