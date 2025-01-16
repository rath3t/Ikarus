// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverlogger.hh
 * \brief Observer implementation for logging non-linear solvers
 */

#pragma once
#include "observer.hh"
#include "observermessages.hh"

#include <ikarus/utils/broadcaster/listener.hh>

namespace Ikarus {
/**
 * \brief Implementation of an observer for logging non-linear solvers.
 * \ingroup observer
 * This class inherits from the IObserver class and provides specific
 * implementations for updating based on NonLinearSolverMessages.
 */
class NonLinearSolverLogger : public Listener
{
public:
  template <typename NLS>
  NonLinearSolverLogger(NLS& nls) {
    // this->subscribe(nls, [&](NonLinearSolverMessages message) { this->updateImpl(message);
    // });
    // this->subscribe(
    //     nls, [&](NonLinearSolverMessages message, double val) { this->updateImpl(message, val); });
    // this->subscribe(
    //     nls, [&](NonLinearSolverMessages message, int intVal) { this->updateImpl(message, intVal); });

    this->subscribe<NLS, void(NonLinearSolverMessages)>(
        nls, [&](NonLinearSolverMessages message) { this->updateImpl(message); });
    this->subscribe<NLS, void(NonLinearSolverMessages, double)>(
        nls, [&](NonLinearSolverMessages message, double val) { this->updateImpl(message, val); });
    this->subscribe<NLS, void(NonLinearSolverMessages, int)>(
        nls, [&](NonLinearSolverMessages message, int intVal) { this->updateImpl(message, intVal); });

    // nls->template station<void(NonLinearSolverMessages)>().registerListener(
    //     [&](NonLinearSolverMessages message) { this->updateImpl(message); });
    // nls->template station<void(NonLinearSolverMessages, double)>().registerListener(
    //     [&](NonLinearSolverMessages message, double val) { this->updateImpl(message, val); });
    // nls->template station<void(NonLinearSolverMessages, int)>().registerListener(
    //     [&](NonLinearSolverMessages message, int val) { this->updateImpl(message, val); });
  }

  /**
   * \brief Handles the update when a NonLinearSolverMessages is received.
   *
   * \param message The NonLinearSolverMessages received.
   */
  void updateImpl(NonLinearSolverMessages message);

  /**
   * \brief Handles the update when a NonLinearSolverMessages with a double value is received.
   *
   * \param message The NonLinearSolverMessages received.
   * \param val The double value associated with the message.
   */
  void updateImpl(NonLinearSolverMessages message, double val);

  /**
   * \brief Handles the update when a NonLinearSolverMessages with an integer value is received.
   *
   * \param message The NonLinearSolverMessages received.
   * \param intVal The integer value associated with numberOfIterations.
   */
  void updateImpl(NonLinearSolverMessages message, int intVal);

private:
  int iters_{0};
  double dNorm_{0};
  double rNorm_{0};
  double lambda_{0};
};
} // namespace Ikarus
