// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controlroutinebase.hh
 * \brief Base for all control routines
 */

#pragma once
#include <ikarus/utils/broadcaster/broadcaster.hh>
#include <ikarus/utils/observer/broadcastermessages.hh>

namespace Ikarus {

struct ControlRoutineBase
    : public Broadcasters<void(ControlMessages), void(ControlMessages, const std::string&),
                          void(ControlMessages, int, const std::string&), void(ControlMessages, int, double)>

{
};

} // namespace Ikarus
