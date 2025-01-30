// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/*!
 * \file loadcontrol.inl
 * \brief Implementation of the run function
 * @ingroup  controlroutines
 */

#pragma once

namespace Ikarus {
template <typename NLS>
ControlInformation LoadControl<NLS>::run() {
  using enum ControlMessages;

  ControlInformation info({false});
  auto& nonOp = nonLinearSolver_->nonLinearOperator();
  this->notify(CONTROL_STARTED, static_cast<std::string>(this->name()));
  auto& loadParameter = nonOp.lastParameter();

  loadParameter = 0.0;
  this->notify(STEP_STARTED, 0, stepSize_);
  auto solverInfo = nonLinearSolver_->solve();
  info.solverInfos.push_back(solverInfo);
  info.totalIterations += solverInfo.iterations;
  if (not solverInfo.success)
    return info;
  this->notify(SOLUTION_CHANGED);
  this->notify(STEP_ENDED);

  for (int ls = 0; ls < loadSteps_; ++ls) {
    this->notify(STEP_STARTED, ls, stepSize_);
    loadParameter += stepSize_;
    solverInfo = nonLinearSolver_->solve();
    info.solverInfos.push_back(solverInfo);
    info.totalIterations += solverInfo.iterations;
    if (not solverInfo.success)
      return info;
    this->notify(SOLUTION_CHANGED);
    this->notify(STEP_ENDED);
  }
  this->notify(CONTROL_ENDED, info.totalIterations, static_cast<std::string>(this->name()));
  info.success = true;
  return info;
}
} // namespace Ikarus
