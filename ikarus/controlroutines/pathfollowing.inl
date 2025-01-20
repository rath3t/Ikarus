// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/*!
 * \file pathfollowing.inl
 * \brief Implementation of the run function
 * @ingroup  controlroutines
 */

#pragma once

#include <Eigen/Core>

#include <ikarus/controlroutines/adaptivestepsizing.hh>
#include <ikarus/controlroutines/controlinfos.hh>
#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphsonwithscalarsubsidiaryfunction.hh>
#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/nonlinearoperator.hh>

namespace Ikarus {

template <typename NLS, typename PF, typename ASS>
requires(Impl::checkPathFollowingTemplates<NLS, PF, ASS>())
ControlInformation PathFollowing<NLS, PF, ASS>::run() {
  ControlInformation info;
  auto& nonOp = nonLinearSolver_->nonLinearOperator();
  this->notifyListeners(ControlMessages::CONTROL_STARTED, pathFollowingType_.name());

  info.totalIterations = 0;
  subsidiaryArgs_.reset();
  subsidiaryArgs_.stepSize = stepSize_;
  subsidiaryArgs_.DD.resizeLike(nonOp.firstParameter());
  subsidiaryArgs_.DD.setZero();
  subsidiaryArgs_.dfdDD.resizeLike(nonOp.firstParameter());
  subsidiaryArgs_.dfdDD.setZero();

  /// Initializing solver
  this->notifyListeners(ControlMessages::STEP_STARTED, 0, subsidiaryArgs_.stepSize);
  pathFollowingType_.initialPrediction(nonOp, subsidiaryArgs_);
  auto solverInfo = nonLinearSolver_->solve(pathFollowingType_, subsidiaryArgs_);
  info.solverInfos.push_back(solverInfo);
  info.totalIterations += solverInfo.iterations;
  if (not solverInfo.success)
    return info;
  this->notifyListeners(ControlMessages::SOLUTION_CHANGED);
  this->notifyListeners(ControlMessages::STEP_ENDED);

  /// Calculate predictor for a particular step
  for (int ls = 1; ls < steps_; ++ls) {
    subsidiaryArgs_.currentStep = ls;

    adaptiveStepSizing_(solverInfo, subsidiaryArgs_, nonOp);

    this->notifyListeners(ControlMessages::STEP_STARTED, subsidiaryArgs_.currentStep, subsidiaryArgs_.stepSize);

    pathFollowingType_.intermediatePrediction(nonOp, subsidiaryArgs_);

    solverInfo = nonLinearSolver_->solve(pathFollowingType_, subsidiaryArgs_);

    info.solverInfos.push_back(solverInfo);
    info.totalIterations += solverInfo.iterations;
    if (not solverInfo.success)
      return info;
    this->notifyListeners(ControlMessages::SOLUTION_CHANGED);
    this->notifyListeners(ControlMessages::STEP_ENDED);
  }

  this->notifyListeners(ControlMessages::CONTROL_ENDED, info.totalIterations, pathFollowingType_.name());
  info.success = true;
  return info;
}
} // namespace Ikarus
