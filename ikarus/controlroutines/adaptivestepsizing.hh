// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file adaptivestepsizing.hh
 * \brief Contains the AdaptiveStepSizing namespace with strategies for adaptive step sizing.
 *
 * This file defines two strategies for adaptive step sizing: NoOp and IterationBased.
 * \ingroup  controlroutines
 */

#pragma once

#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/solver/nonlinearsolver/solverinfos.hh>

namespace Ikarus::AdaptiveStepSizing {

/**
 * \brief The NoOp strategy for adaptive step sizing.
 *
 * This strategy does nothing and keeps the step size constant.
 */
struct NoOp
{
  /**
   * \brief Call operator for the NoOp strategy.
   *
   * Does nothing.
   *
   * \param solverInfo Information about the nonlinear solver.
   * \param subsidiaryArgs Subsidiary arguments for adaptive step sizing.
   * \param nonLinearOperator The nonlinear operator.
   * \tparam NLO The nonlinear operator type.
   */
  template <typename NLO>
  void operator()(const NonLinearSolverInformation& solverInfo, SubsidiaryArgs& subsidiaryArgs,
                  const NLO& nonLinearOperator) {}

  /**
   * \brief Get the target iterations.
   * \return Always returns 0, indicating no specific target.
   */
  [[nodiscard]] int targetIterations() const { return 0; }

  /**
   * \brief Set the target iterations.
   * \param targetIterations The target iterations (ignored).
   */
  void setTargetIterations([[maybe_unused]] int targetIterations) {}
};

/**
 * \brief The IterationBased strategy for adaptive step sizing.
 *
 * This strategy adjusts the step size based on the ratio of target iterations to previous iterations.
 *
 * \details
 * The step size is changed according to the number of iterations of the previous step
 * \f[ \hat{s}_{k+1} = \sqrt{\frac{\hat{i}}{i_k}}\hat{s}_k,  \f] where \f$\hat{s}_{k+1} \f$ is the new stesize and
 * \f$\hat{i}\f$ is the desired number of iterations and \f$i_k\f$ is the number of iterations of the previous step.
 */
struct IterationBased
{
  /**
   * \brief Call operator of the IterationBased strategy.
   *
   * Adjusts the step size based on the ratio of target iterations to previous iterations.
   *
   * \param solverInfo Information about the nonlinear solver.
   * \param subsidiaryArgs Subsidiary arguments for adaptive step sizing.
   * \param nonLinearOperator The nonlinear operator.
   * \tparam NLO The nonlinear operator.
   */
  template <typename NLO>
  void operator()(const NonLinearSolverInformation& solverInfo, SubsidiaryArgs& subsidiaryArgs,
                  const NLO& nonLinearOperator) {
    if (subsidiaryArgs.currentStep == 0)
      return;
    if (targetIterations_ == 0)
      DUNE_THROW(Dune::InvalidStateException,
                 "TargetIterations should not be equal to 0. Try calling "
                 "setTargetIterations(int) first.");
    const auto previousIterations = solverInfo.iterations;
    if (previousIterations > targetIterations_)
      subsidiaryArgs.stepSize =
          sqrt(static_cast<double>(targetIterations_) / previousIterations) * subsidiaryArgs.stepSize;
  }

  /**
   * \brief Get the target iterations.
   * \return The number of target iterations.
   */
  [[nodiscard]] int targetIterations() const { return targetIterations_; }

  /**
   * \brief Set the target iterations.
   * \param targetIterations The number of target iterations.
   */
  void setTargetIterations(int targetIterations) {
    if (targetIterations == 0)
      DUNE_THROW(Dune::InvalidStateException, "TargetIterations should not be equal to 0.");
    targetIterations_ = targetIterations;
  }

private:
  int targetIterations_{0};
};
} // namespace Ikarus::AdaptiveStepSizing
