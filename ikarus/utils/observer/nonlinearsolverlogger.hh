// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverlogger.hh
 * \brief Observer implementation for logging non-linear solvers
 */

#pragma once
#include "observer.hh"
#include "observermessages.hh"

namespace Ikarus {
/**
 * \brief Implementation of an observer for logging non-linear solvers.
 * \ingroup observer
 * This class inherits from the IObserver class and provides specific
 * implementations for updating based on NonLinearSolverMessages.
 */
class NonLinearSolverLogger : public IObserver<NonLinearSolverMessages, Eigen::VectorXd, Eigen::VectorXd>
{
public:
  /**
   * \brief Handles the update when a NonLinearSolverMessages is received.
   *
   * \param message The NonLinearSolverMessages received.
   */
  void updateImpl(NonLinearSolverMessages message) final;

  /**
   * \brief Handles the update when a NonLinearSolverMessages with a double value is received.
   *
   * \param message The NonLinearSolverMessages received.
   * \param val The double value associated with the message.
   */
  void updateImpl(NonLinearSolverMessages message, double val) final;

  /**
   * \brief Handles the update when a NonLinearSolverMessages with an integer value is received.
   *
   * \param message The NonLinearSolverMessages received.
   * \param intVal The integer value associated with numberOfIterations.
   */
  void updateImpl(NonLinearSolverMessages message, int intVal) final;

private:
  int iters_{0};
  double dNorm_{0};
  double rNorm_{0};
  double lambda_{0};
};
} // namespace Ikarus
