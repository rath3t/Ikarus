// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file convergencecriteria.hh
 * \brief A collection of convergence criteria used by the nonlinear solvers.
 */

#include <dune/common/float_cmp.hh>

#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/utils/defaultfunctions.hh>
#include <ikarus/utils/linearalgebrahelper.hh>
#include <ikarus/utils/makeenum.hh>

#pragma once

namespace Ikarus {

/**
 * \enum PreConditioner
 * \brief Enumeration of available preconditioners for the trust region solver.
 */
enum class PreConditioner
{
  IncompleteCholesky,
  IdentityPreconditioner,
  DiagonalPreconditioner
};

MAKE_ENUM(ConvergenceCriterion, ResiduumNorm, CorrectionNorm, MaximumOfCorrectionVector);

template <typename NLO, ConvergenceCriterion CC = ConvergenceCriterion::ResiduumNorm,
          typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
class NewtonRaphson;

template <typename NLO, ConvergenceCriterion CC = ConvergenceCriterion::ResiduumNorm,
          typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
class NewtonRaphsonWithSubsidiaryFunction;

template <typename NLO, PreConditioner preConditioner = PreConditioner::IncompleteCholesky,
          typename UF = utils::UpdateDefault>
class TrustRegion;

template <ConvergenceCriterion CC>
requires(CC != ConvergenceCriterion::BEGIN and CC != ConvergenceCriterion::END)
struct NRConvergenceCriteria
{
  template <typename NLO, typename NSS, typename CV>
  bool operator()(const NLO& nonLinearOperator, const NSS& solverSettings, const CV& correctionVector) {
    if constexpr (CC == ConvergenceCriterion::ResiduumNorm) {
      const auto& rx = nonLinearOperator.value();
      auto rNorm     = norm(rx);
      return Dune::FloatCmp::le(static_cast<double>(rNorm), solverSettings.tol);
    } else if constexpr (CC == ConvergenceCriterion::CorrectionNorm) {
      auto dNorm = norm(correctionVector);
      return Dune::FloatCmp::le(static_cast<double>(dNorm), solverSettings.tol);
    } else if constexpr (CC == ConvergenceCriterion::MaximumOfCorrectionVector) {
      auto maxCorr = max(correctionVector);
      return Dune::FloatCmp::le(static_cast<double>(maxCorr), solverSettings.tol);
    } else {
      return false; // TODO Change to Dune::AlwaysFalse
    }
  }
};

template <ConvergenceCriterion CC>
requires(CC != ConvergenceCriterion::BEGIN and CC != ConvergenceCriterion::END)
struct NRWSFConvergenceCriteria
{
  template <typename NLO, typename NSS, typename DD>
  bool operator()(const NLO& nonLinearOperator, const NSS& solverSettings, const DD& deltaD, double deltaLambda,
                  const SubsidiaryArgs& subsidiaryArgs) {
    if constexpr (CC == ConvergenceCriterion::ResiduumNorm) {
      const auto& rx = nonLinearOperator.value();
      Eigen::VectorXd totalResidual(rx.size() + 1);
      totalResidual << rx, subsidiaryArgs.f;
      auto rNorm = norm(totalResidual);
      return Dune::FloatCmp::le(static_cast<double>(rNorm), solverSettings.tol);
    } else if constexpr (CC == ConvergenceCriterion::CorrectionNorm) {
      Eigen::VectorXd correctionVector(deltaD.size() + 1);
      correctionVector << deltaD, deltaLambda;
      auto dNorm = norm(correctionVector);
      return Dune::FloatCmp::le(static_cast<double>(dNorm), solverSettings.tol);
    } else if constexpr (CC == ConvergenceCriterion::MaximumOfCorrectionVector) {
      Eigen::VectorXd correctionVector(deltaD.size() + 1);
      correctionVector << deltaD, deltaLambda;
      auto maxCorr = max(correctionVector);
      return Dune::FloatCmp::le(static_cast<double>(maxCorr), solverSettings.tol);
    } else {
      return false; // TODO Change to Dune::AlwaysFalse
    }
  }
};

template <ConvergenceCriterion CC>
requires(CC != ConvergenceCriterion::BEGIN and CC != ConvergenceCriterion::END)
struct TRConvergenceCriteria
{
  template <typename NLO, typename NSS, typename CV>
  bool operator()(const NLO& nonLinearOperator, const NSS& solverSettings, const CV& correctionVector) {
    if constexpr (CC == ConvergenceCriterion::ResiduumNorm) {
      return false;
    } else if constexpr (CC == ConvergenceCriterion::CorrectionNorm) {
      return false;
    } else if constexpr (CC == ConvergenceCriterion::MaximumOfCorrectionVector) {
      return false;
    } else {
      return false; // TODO Change to Dune::AlwaysFalse
    }
  }
};

template <typename NLS, ConvergenceCriterion CC>
struct ConvergenceCriteria
{
  using NR    = NewtonRaphson<typename NLS::NonLinearOperator, NLS::criteraType, typename NLS::LinearSolverType,
                              typename NLS::UpdateFunctionType>;
  using NRWSF = NewtonRaphsonWithSubsidiaryFunction<typename NLS::NonLinearOperator, NLS::criteraType,
                                                    typename NLS::LinearSolverType, typename NLS::UpdateFunctionType>;
  using TR    = TrustRegion<typename NLS::NonLinearOperator>;

  using Criteria = std::conditional_t<
      std::is_same_v<NLS, NR>, NRConvergenceCriteria<CC>,
      std::conditional_t<std::is_same_v<NLS, NRWSF>, NRWSFConvergenceCriteria<CC>, TRConvergenceCriteria<CC>>>;

  template <typename... Args>
  bool operator()(Args&&... args) {
    return criteria(std::forward<Args>(args)...);
  }

private:
  Criteria criteria{};
};

} // namespace Ikarus
