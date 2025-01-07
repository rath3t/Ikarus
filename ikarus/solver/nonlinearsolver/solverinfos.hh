// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file solverinfos.hh
 * \brief Implementation of the Newton-Raphson method for solving nonlinear equations.
 */

#pragma once

namespace Ikarus {
/**
 * \brief Information about the result of a non-linear solver.
 *
 * This structure holds information about the success of a non-linear solver,
 * including the residual norm, correction norm, and the number of iterations performed.
 */
struct NonLinearSolverInformation
{
  /**
   * \brief Conversion operator to check the success of the solver.
   *
   * \return `true` if the solver was successful, `false` otherwise.
   */
  explicit operator bool() const { return success; }
  bool success{false};
  double residualNorm{std::numeric_limits<double>::infinity()};
  double correctionNorm{std::numeric_limits<double>::infinity()};
  int iterations{-1};
};
} // namespace Ikarus
