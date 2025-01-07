// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file observermessages.hh
 * \brief Enums for observer messages
 */

#pragma once

namespace Ikarus {
/**
 * \brief Enum class defining control-routine-related messages.
 * \ingroup observer
 */
enum class ControlMessages
{
  BEGIN,
  CONTROL_STARTED,
  CONTROL_ENDED,
  STEP_STARTED,
  STEP_ENDED,
  SOLUTION_CHANGED,
  END
};

/**
 * \brief Enum class defining non-linear solver-related messages.
 * \ingroup observer
 */
enum class NonLinearSolverMessages
{
  BEGIN,
  INIT,
  ITERATION_STARTED,
  ITERATION_ENDED,
  RESIDUALNORM_UPDATED,
  CORRECTIONNORM_UPDATED,
  SOLUTION_CHANGED,
  FINISHED_SUCESSFULLY,
  END
};
} // namespace Ikarus
