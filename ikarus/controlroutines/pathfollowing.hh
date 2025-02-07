// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file pathfollowing.hh
 * \brief Defines the PathFollowing class.
 */

#pragma once

#include <memory>

#include <ikarus/controlroutines/adaptivestepsizing.hh>
#include <ikarus/controlroutines/controlinfos.hh>
#include <ikarus/controlroutines/controlroutinebase.hh>
#include <ikarus/controlroutines/controlroutinefactory.hh>
#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphsonwithscalarsubsidiaryfunction.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverfactory.hh>
#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/nonlinopfactory.hh>

namespace Ikarus {

namespace Impl {

  /**
   * \brief Checks template requirements for path-following control routines.
   *
   * This consteval function checks whether the provided template parameters meet the requirements
   * for path-following control routines.
   *
   * \tparam NLS Type of the nonlinear solver.
   * \tparam PF Type of the path-following strategy.
   * \tparam ASS Type of the adaptive step sizing strategy.
   *
   * \return True if the templates meet the requirements, false otherwise.
   */
  template <typename NLS, typename PF = ArcLength, typename ASS>
  consteval bool checkPathFollowingTemplates() {
    return Concepts::PathFollowingStrategy<PF, typename NLS::NonLinearOperator, SubsidiaryArgs> and
           Concepts::AdaptiveStepSizingStrategy<ASS, NonLinearSolverInformation, SubsidiaryArgs,
                                                std::remove_cvref_t<typename NLS::NonLinearOperator>> and
           Concepts::NonLinearSolverCheckForPathFollowing<NLS>;
  }

} // namespace Impl

template <typename NLS, typename PF, typename ASS>
requires(Impl::checkPathFollowingTemplates<NLS, PF, ASS>())
class PathFollowing;

/**
 * \struct PathFollowingConfig
 * \brief Config for the Path-Following control routine
 *
 * \tparam PF_ the type of PathFollowing that is used (defaults to ArcLength)
 * \tparam ASS_ the type of AdaptiveStepSizing that is used (defaults to NoOp)
 */
template <typename PF_ = ArcLength, typename ASS_ = AdaptiveStepSizing::NoOp>
struct PathFollowingConfig
{
  using PF  = PF_;
  using ASS = ASS_;

  int steps{};
  double stepSize{};
  PF pathFollowingFunction{};
  ASS adaptiveStepSizingFunction{};
};

#ifndef DOXYGEN
PathFollowingConfig(int, double) -> PathFollowingConfig<>;

template <typename PF>
PathFollowingConfig(int, double, PF) -> PathFollowingConfig<PF>;

template <typename PF, typename ASS>
PathFollowingConfig(int, double, PF, ASS) -> PathFollowingConfig<PF, ASS>;

template <typename ASS>
PathFollowingConfig(int, double, ArcLength, ASS) -> PathFollowingConfig<ArcLength, ASS>;
#endif

/**
 * \brief Function to create a path following instance
 *
 * \tparam NLS Type of the nonlinear solver
 * \tparam PFConfig  the provided config for the path following
 * \param config the provided config for the path following
 * \param nonlinearSolver the provided nonlinearsolver
 * \return PathFollowing
 */
template <typename NLS, typename PFConfig>
requires traits::isSpecialization<PathFollowingConfig, std::remove_cvref_t<PFConfig>>::value
auto createControlRoutine(PFConfig&& config, NLS&& nonlinearSolver) {
  return PathFollowing<typename std::remove_cvref_t<NLS>::element_type, typename PFConfig::PF, typename PFConfig::ASS>(
      std::forward<NLS>(nonlinearSolver), config.steps, config.stepSize, config.pathFollowingFunction,
      config.adaptiveStepSizingFunction);
}

/**
 * \class PathFollowing
 * \brief The PathFollowing control routine for path-following analysis.
 *
 * This class represents the PathFollowing control routine, which utilizes a nonlinear solver,
 * such as Newton's method with scalar subsidiary function, which has to be fulfilled for solving the system along a
 * predefined path.
 *
 * \details Consider a non-linear system of equations
 *        \f[\mathbf{R}: \require{cases}\begin{cases}\mathbb{R}^n \times \mathbb{R} &\rightarrow \mathbb{R}^n
 * \\ (\mathbf{D},\lambda)
 * &\mapsto \mathbf{R}(\mathbf{D},\lambda) \end{cases}.\f]
 *
 * Then in each step \f$k+1\f$ of the path following algorithm, the following problem is solved
 * \f[ \begin{align}
 * \mathbf{R}(\mathbf{D}_k+ \mathrm{D}\mathbf{D}, \lambda_k+ \mathrm{D} \lambda) &= \mathbf{0} \\
 * f(\mathrm{D}\mathbf{D}, \mathrm{D} \lambda) &= 0 \end{align} \f]
 *
 * where \f$\mathrm{D}\mathbf{D}\f$ is the increment of the solution vector and \f$\mathrm{D} \lambda\f$ is the load
 * factor increment. The subsidiary function \f$f\f$ is provided by the user and needs to fulfill the concept
 * Concepts::PathFollowingStrategy. This subsidiary function makes the given problem well-posed.
 *
 *  Currently the following subsidiary functions are implemented \ref LoadControlSubsidiaryFunction, ArcLength and
 * DisplacementControl
 * \ingroup controlroutines
 * \tparam NLS Type of the nonlinear solver used in the control routine.
 * \tparam PF Type of the path-following strategy.
 * \tparam ASS Type of the adaptive step sizing strategy.
 */
template <typename NLS, typename PF = ArcLength, typename ASS = AdaptiveStepSizing::NoOp>
requires(Impl::checkPathFollowingTemplates<NLS, PF, ASS>())
class PathFollowing : public ControlRoutineBase<typename NLS::NonLinearOperator>

{
public:
  /** \brief The name of the PathFollowing method. */
  constexpr auto name() const { return std::string("Path following with " + pathFollowingType_.name()); }

  /**
   * \brief Constructor for PathFollowing.
   * \param nls The non linear solver.
   * \param steps Number of steps in the control routine.
   * \param stepSize Size of each step.
   * \param pathFollowingType Type of the path-following function.
   * \param adaptiveStepSizing Type of the adaptive step sizing strategy.
   */
  PathFollowing(const std::shared_ptr<NLS>& nls, int steps, double stepSize, PF pathFollowingType = ArcLength{},
                ASS adaptiveStepSizing = {})
      : nonLinearSolver_{nls},
        steps_{steps},
        stepSize_{stepSize},
        pathFollowingType_{pathFollowingType},
        adaptiveStepSizing_{adaptiveStepSizing} {}

  /**
   * \brief Executes the PathFollowing routine.
   *
   * \return ControlInformation structure containing information about the control results.
   */
  ControlInformation run();

  /* \brief returns the nonlinear solver */
  NLS& nonlinearSolver() { return *nonLinearSolver_; }

private:
  std::shared_ptr<NLS> nonLinearSolver_;

  int steps_;
  double stepSize_;
  PF pathFollowingType_;
  ASS adaptiveStepSizing_;
  SubsidiaryArgs subsidiaryArgs_;
};

} // namespace Ikarus

#include <ikarus/controlroutines/pathfollowing.inl>
