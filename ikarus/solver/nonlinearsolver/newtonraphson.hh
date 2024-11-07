// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file newtonraphson.hh
 * \brief Implementation of the Newton-Raphson method for solving nonlinear equations.
 */

#pragma once

#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/solver/nonlinearsolver/helperfunctions.hh>
#include <ikarus/solver/nonlinearsolver/solverinfos.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/defaultfunctions.hh>
#include <ikarus/utils/linearalgebrahelper.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

template <typename NLO, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
class NewtonRaphson;

struct NRSettings
{
  double tol{1e-8};
  int maxIter{20};
  int minIter{0};
};

/**
 * \struct NewtonRaphsonConfig
 * \brief Config for the Newton-Raphson solver.
 */
template <typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
struct NewtonRaphsonConfig
{
  using LinearSolver   = LS;
  using UpdateFunction = UF;
  NRSettings parameters;
  LS linearSolver;
  UF updateFunction;

  template <typename UF2>
  auto rebindUpdateFunction(UF2&& updateFunction) const {
    NewtonRaphsonConfig<LS, UF2> settings{
        .parameters = parameters, .linearSolver = linearSolver, .updateFunction = std::forward<UF2>(updateFunction)};
    return settings;
  }

  template <typename NLO>
  using Solver = NewtonRaphson<NLO, LS, UF>;
};

/**
 * \brief Function to create a NewtonRaphson solver instance.
 * \tparam NLO Type of the nonlinear operator to solve.
 * \tparam NRConfig Type of the nonlinear solver config.
 * \param config Config for the solver.
 * \param nonLinearOperator Nonlinear operator to solve.
 * \return Shared pointer to the NewtonRaphson solver instance.
 */
template <typename NLO, typename NRConfig>
requires traits::isSpecialization<NewtonRaphsonConfig, std::remove_cvref_t<NRConfig>>::value
auto createNonlinearSolver(NRConfig&& config, NLO&& nonLinearOperator) {
  using LS           = std::remove_cvref_t<NRConfig>::LinearSolver;
  using UF           = std::remove_cvref_t<NRConfig>::UpdateFunction;
  auto solverFactory = []<class NLO2, class LS2, class UF2>(NLO2&& nlo2, LS2&& ls, UF2&& uf) {
    return std::make_shared<
        NewtonRaphson<std::remove_cvref_t<NLO2>, std::remove_cvref_t<LS2>, std::remove_cvref_t<UF2>>>(
        nlo2, std::forward<LS2>(ls), std::forward<UF2>(uf));
  };

  if constexpr (std::remove_cvref_t<NLO>::numberOfFunctions == 3) {
    auto solver =
        solverFactory(nonLinearOperator.template subOperator<1, 2>(), std::forward<NRConfig>(config).linearSolver,
                      std::forward<NRConfig>(config).updateFunction);
    solver->setup(config.parameters);
    return solver;
  } else {
    static_assert(std::remove_cvref_t<NLO>::numberOfFunctions > 1,
                  "The number of derivatives in the nonlinear operator have to be more than 1");
    auto solver = solverFactory(nonLinearOperator, std::forward<NRConfig>(config).linearSolver,
                                std::forward<NRConfig>(config).updateFunction);
    ;

    solver->setup(std::forward<NRConfig>(config).parameters);
    return solver;
  }
}

/**
 * \class NewtonRaphson
 * \brief Implementation of the Newton-Raphson method for solving nonlinear equations.
 * \tparam NLO Type of the nonlinear operator to solve.
 * \tparam LS Type of the linear solver used internally (default is SolverDefault).
 * \tparam UF Type of the update function (default is UpdateDefault).
 * \relates makeNewtonRaphson
 * \ingroup solvers
 */
template <typename NLO, typename LS, typename UF>
class NewtonRaphson : public IObservable<NonLinearSolverMessages>
{
public:
  using Settings = NRSettings;
  ///< Compile-time boolean indicating if the linear solver satisfies the non-linear solver concept
  static constexpr bool isLinearSolver =
      Ikarus::Concepts::LinearSolverCheck<LinearSolver, typename NLO::DerivativeType, typename NLO::ValueType>;

  ///< Type representing the parameter vector of the nonlinear operator.
  using ValueType = typename NLO::template ParameterValue<0>;

  using UpdateFunction    = UF;  ///< Type representing the update function.
  using NonLinearOperator = NLO; ///< Type of the non-linear operator

  /**
   * \brief Constructor for NewtonRaphson.
   * \param nonLinearOperator Nonlinear operator to solve.
   * \param linearSolver Linear solver used internally (default is SolverDefault).
   * \param updateFunction Update function (default is UpdateDefault).
   */
  template <typename LS2 = LS, typename UF2 = UF>
  explicit NewtonRaphson(const NonLinearOperator& nonLinearOperator, LS2&& linearSolver = {}, UF2&& updateFunction = {})
      : nonLinearOperator_{nonLinearOperator},
        linearSolver_{std::forward<LS2>(linearSolver)},
        updateFunction_{std::forward<UF2>(updateFunction)} {
    if constexpr (std::is_same_v<typename NonLinearOperator::ValueType, Eigen::VectorXd>)
      correction_.setZero(this->nonLinearOperator().value().size());
  }

  /**
   * \brief Set up the solver with the given settings.
   * \param settings Newton-Raphson settings.
   */
  void setup(const Settings& settings) {
    if (settings.minIter > settings.maxIter)
      DUNE_THROW(Dune::InvalidStateException,
                 "Minimum number of iterations cannot be greater than maximum number of iterations");
    settings_ = settings;
  }

#ifndef DOXYGEN
  struct NoPredictor
  {
  };
#endif
  /**
   * \brief Solve the nonlinear system.
   * \param dxPredictor Predictor for the solution increment (default is NoPredictor).
   * \return Information about the solution process.
   */
  template <typename SolutionType = NoPredictor>
  requires std::is_same_v<SolutionType, NoPredictor> ||
           std::is_convertible_v<SolutionType, std::remove_cvref_t<typename NonLinearOperator::ValueType>>
  [[nodiscard(
      "The solve method returns information of the solution process. You should store this information and check if "
      "it was successful")]] Ikarus::NonLinearSolverInformation
  solve(const SolutionType& dxPredictor = NoPredictor{}) {
    this->notify(NonLinearSolverMessages::INIT);
    Ikarus::NonLinearSolverInformation solverInformation;
    solverInformation.success = true;
    auto& x                   = nonLinearOperator().firstParameter();
    if constexpr (not std::is_same_v<SolutionType, NoPredictor>)
      updateFunction_(x, dxPredictor);
    nonLinearOperator().updateAll();
    const auto& rx = nonLinearOperator().value();
    const auto& Ax = nonLinearOperator().derivative();
    auto rNorm     = norm(rx);
    decltype(rNorm) dNorm;
    int iter{0};
    if constexpr (isLinearSolver)
      linearSolver_.analyzePattern(Ax);
    while ((rNorm > settings_.tol && iter < settings_.maxIter) or iter < settings_.minIter) {
      this->notify(NonLinearSolverMessages::ITERATION_STARTED);
      if constexpr (isLinearSolver) {
        linearSolver_.factorize(Ax);
        linearSolver_.solve(correction_, -rx);
        dNorm = correction_.norm();
      } else {
        correction_ = -linearSolver_(rx, Ax);
        dNorm       = norm(correction_);
      }
      updateStates(nonLinearOperator().assembler(), correction_);
      updateFunction_(x, correction_);
      this->notify(NonLinearSolverMessages::CORRECTIONNORM_UPDATED, static_cast<double>(dNorm));
      this->notify(NonLinearSolverMessages::SOLUTION_CHANGED);
      nonLinearOperator().updateAll();
      rNorm = norm(rx);
      this->notify(NonLinearSolverMessages::RESIDUALNORM_UPDATED, static_cast<double>(rNorm));
      this->notify(NonLinearSolverMessages::ITERATION_ENDED);
      ++iter;
    }
    if (iter == settings_.maxIter)
      solverInformation.success = false;
    solverInformation.iterations     = iter;
    solverInformation.residualNorm   = static_cast<double>(rNorm);
    solverInformation.correctionNorm = static_cast<double>(dNorm);
    if (solverInformation.success)
      this->notify(NonLinearSolverMessages::FINISHED_SUCESSFULLY, iter);
    return solverInformation;
  }

  /**
   * \brief Access the nonlinear operator.
   * \return Reference to the nonlinear operator.
   */
  auto& nonLinearOperator() { return nonLinearOperator_; }

private:
  NonLinearOperator nonLinearOperator_;
  typename NonLinearOperator::ValueType correction_;
  LS linearSolver_;
  UpdateFunction updateFunction_;
  Settings settings_;
};

/**
 * \brief Function to create a NewtonRaphson solver instance.
 * \tparam NLO Type of the nonlinear operator to solve.
 * \tparam LS Type of the linear solver used internally (default is SolverDefault).
 * \tparam UF Type of the update function (default is UpdateDefault).
 * \param nonLinearOperator Nonlinear operator to solve.
 * \param linearSolver Linear solver used internally (default is SolverDefault).
 * \param updateFunction Update function (default is UpdateDefault).
 * \return Shared pointer to the NewtonRaphson solver instance.
 */
template <typename NLO, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
auto makeNewtonRaphson(const NLO& nonLinearOperator, LS&& linearSolver = {}, UF&& updateFunction = {}) {
  return std::make_shared<NewtonRaphson<NLO, LS, UF>>(nonLinearOperator, std::forward<LS>(linearSolver),
                                                      std::move(updateFunction));
}

template <typename NLO, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
NewtonRaphson(const NLO& nonLinearOperator, LS&& linearSolver = {},
              UF&& updateFunction = {}) -> NewtonRaphson<NLO, std::remove_cvref_t<LS>, std::remove_cvref_t<UF>>;

} // namespace Ikarus
