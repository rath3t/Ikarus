// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "nonlinearsolverlogger.hh"

#include "spdlog/spdlog.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wswitch-enum"

namespace Ikarus {

void NonLinearSolverLogger::updateImpl(NonLinearSolverMessages message) {
  switch (message) {
    case NonLinearSolverMessages::ITERATION_STARTED:
      break;
    case NonLinearSolverMessages::INIT:
      iters_ = 1;
      rNorm_ = 0.0;
      dNorm_ = 0.0;
      spdlog::info("Non-linear solver started:");
      spdlog::info("{:<11} {:<20} {:<20} {:<20}", "Ite", "normR", "normD", "lambda");
      spdlog::info("-------------------------------------------------------------------------------");
      break;
    case NonLinearSolverMessages::ITERATION_ENDED:
      spdlog::info("{} {:<10d} {:<20.2e} {:<20.2e} {:<20.2e}", "", iters_, rNorm_, dNorm_, lambda_);
      ++iters_;
      break;
    default:
      break;
  }
}

void NonLinearSolverLogger::updateImpl(NonLinearSolverMessages message, double val) {
  switch (message) {
    case NonLinearSolverMessages::RESIDUALNORM_UPDATED:
      rNorm_ = val;
      break;
    case NonLinearSolverMessages::SOLUTION_CHANGED:
      lambda_ = val;
      break;
    case NonLinearSolverMessages::CORRECTIONNORM_UPDATED:
      dNorm_ = val;
      break;
    default:
      break;
  }
}

void NonLinearSolverLogger::updateImpl(NonLinearSolverMessages message, int numberOfIterations) {
  switch (message) {
    case NonLinearSolverMessages::FINISHED_SUCESSFULLY:
      spdlog::info("Number of iterations: {}", numberOfIterations);
      break;
    default:
      break;
  }
}
} // namespace Ikarus

#pragma GCC diagnostic pop
