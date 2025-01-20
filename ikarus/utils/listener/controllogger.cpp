// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "controllogger.hh"

#include "spdlog/spdlog.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wswitch-enum"

namespace Ikarus {

void ControlLogger::updateImpl(ControlMessages message) {
  switch (message) {
    case ControlMessages::STEP_ENDED:
      spdlog::info("===============================================================================");
      break;
    default:
      break; // default: do nothing when notified
  }
}

void ControlLogger::updateImpl(ControlMessages message, const std::string& pathFollowingName) {
  switch (message) {
    case ControlMessages::CONTROL_STARTED:
      start_ = std::chrono::high_resolution_clock::now();
      spdlog::info("===============================================================================");
      spdlog::info("Started path following with: {}", pathFollowingName);
      spdlog::info("===============================================================================");
      break;
    default:
      break;
  }
}

void ControlLogger::updateImpl(ControlMessages message, int stepNumber, double stepSize) {
  switch (message) {
    case ControlMessages::STEP_STARTED:
      spdlog::info("Load step: {:>4} {:>49} {:<.2e}", stepNumber, "Step size = ", stepSize);
      spdlog::info("-------------------------------------------------------------------------------");
      break;
    default:
      break;
  }
}

void ControlLogger::updateImpl(ControlMessages message, int totalIterations, const std::string& pathFollowingName) {
  switch (message) {
    case ControlMessages::CONTROL_ENDED:
      stop_     = std::chrono::high_resolution_clock::now();
      duration_ = duration_cast<std::chrono::milliseconds>(stop_ - start_);
      spdlog::info("End of path following with {} control", pathFollowingName);
      spdlog::info("Total number of iterations: {:3d}", totalIterations);
      spdlog::info("Elapsed time: {} ms", duration_.count());
      break;
    default:
      break;
  }
}
} // namespace Ikarus
#pragma GCC diagnostic pop
