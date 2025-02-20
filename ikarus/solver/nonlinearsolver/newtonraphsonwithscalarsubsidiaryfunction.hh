// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file newtonraphson.hh
 * \brief Implementation of the Newton-Raphson method for solving nonlinear equations.
 */

#pragma once

#include <iosfwd>
#include <type_traits>
#include <utility>

#include "ikarus/assembler/dirichletbcenforcement.hh"
#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/solver/nonlinearsolver/solverinfos.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

template <typename NLO, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
class NewtonRaphsonWithSubsidiaryFunction;

struct NewtonRaphsonWithSubsidiaryFunctionSettings
{
  double tol{1e-8};
  int maxIter{20};
};

/**
 * \struct NewtonRaphsonWithSubsidiaryFunctionConfig
 * \brief Settings for the Newton-Raphson solver with subsidiary function.
 */
template <typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
struct NewtonRaphsonWithSubsidiaryFunctionConfig
{
  using LinearSolver   = LS;
  using UpdateFunction = UF;
  NewtonRaphsonWithSubsidiaryFunctionSettings parameters;
  LS linearSolver;
  UF updateFunction;

  template <typename UF2>
  auto rebindUpdateFunction(UF2&& updateFunction) const {
    NewtonRaphsonWithSubsidiaryFunctionConfig<LS, UF2> settings{
        .parameters = parameters, .linearSolver = linearSolver, .updateFunction = std::forward<UF2>(updateFunction)};
    return settings;
  }

  template <typename NLO>
  using Solver = NewtonRaphsonWithSubsidiaryFunction<NLO, LS, UF>;
};

/**
 * \brief Function to create a NewtonRaphson with subsidiary function solver instance.
 * \tparam NLO Type of the nonlinear operator to solve.
 * \tparam NRConfig Type of the nonlinear solver config.
 * \param config Config for the solver.
 * \param nonLinearOperator Nonlinear operator to solve.
 * \return Shared pointer to the NewtonRaphson solver instance.
 */
template <typename NLO, typename NRConfig>
requires traits::isSpecialization<NewtonRaphsonWithSubsidiaryFunctionConfig, std::remove_cvref_t<NRConfig>>::value
auto createNonlinearSolver(NRConfig&& config, NLO&& nonLinearOperator) {
  using LS           = std::remove_cvref_t<NRConfig>::LinearSolver;
  using UF           = std::remove_cvref_t<NRConfig>::UpdateFunction;
  auto solverFactory = []<class NLO2, class LS2, class UF2>(NLO2&& nlo2, LS2&& ls, UF2&& uf) {
    return std::make_shared<NewtonRaphsonWithSubsidiaryFunction<std::remove_cvref_t<NLO2>, std::remove_cvref_t<LS2>,
                                                                std::remove_cvref_t<UF2>>>(nlo2, std::forward<LS2>(ls),
                                                                                           std::forward<UF2>(uf));
  };

  if constexpr (std::remove_cvref_t<NLO>::nDerivatives == 3) {
    auto solver =
        solverFactory(derivative(nonLinearOperator), std::forward<NRConfig>(config).linearSolver,
                      std::forward<NRConfig>(config).updateFunction);
    solver->setup(config.parameters);
    return solver;
  } else {
    static_assert(std::remove_cvref_t<NLO>::nDerivatives > 1,
                  "The number of derivatives in the nonlinear operator have to be more than 1");
    auto solver = solverFactory(nonLinearOperator, std::forward<NRConfig>(config).linearSolver,
                                std::forward<NRConfig>(config).updateFunction);
    ;

    solver->setup(std::forward<NRConfig>(config).parameters);
    return solver;
  }
}

/**
 * \brief Newton-Raphson solver with subsidiary function.
 *
 * This class provides a Newton-Raphson solver for solving nonlinear systems with a subsidiary function.
 * It uses a linear solver to handle the linear system arising in each iteration.
 *
 * \tparam NLO Type of the nonlinear operator.
 * \tparam LS Type of the linear solver used internally (default is SolverDefault).
 * \tparam UF Type of the update function (default is UpdateDefault).
 */
template <typename NLO, typename LS, typename UF>
class NewtonRaphsonWithSubsidiaryFunction : public IObservable<NonLinearSolverMessages>
{
public:
  using Settings = NewtonRaphsonWithSubsidiaryFunctionSettings;
   using SignatureTraits= typename NLO::Traits;
     ///< Type representing the parameter vector of the nonlinear operator.
  using Domain = typename SignatureTraits::Domain;
  using CorrectionType = std::remove_cvref_t<typename SignatureTraits::template Range<0>>;        ///< Type of the correction of x += deltaX.
  ///< Compile-time boolean indicating if the linear solver satisfies the non-linear solver concept
  static constexpr bool isLinearSolver =
      Ikarus::Concepts::LinearSolverCheck<LS, CorrectionType, Domain>;


  ///< Type representing the update function.
  using UpdateFunctionType = UF;
  using NonLinearOperator  = NLO; ///< Type of the non-linear operator

  using JacobianType = typename SignatureTraits::template Range<1>;

  /**
   * \brief Constructor for NewtonRaphsonWithSubsidiaryFunction.
   * \param residual residual to solve.
   * \param linearSolver Linear solver used internally (default is SolverDefault).
   * \param updateFunction Update function (default is UpdateDefault).
   */
  template <typename LS2 = LS, typename UF2 = UF>
  explicit NewtonRaphsonWithSubsidiaryFunction(const NonLinearOperator& residual, LS2&& linearSolver = {}, UF2&& updateFunction = {})
      : residualFunction_{residual}, jacobianFunction_{derivative(residualFunction_)},
        linearSolver_{std::forward<LS2>(linearSolver)},
        updateFunction_{std::forward<UF2>(updateFunction)} {}

  /**
   * \brief Setup the Newton-Raphson solver with subsidiary function.
   *
   * \param p_settings Settings for the solver.
   */
  void setup(const Settings& settings) { settings_ = settings; }

  /**
   * \brief Solve the nonlinear system using the Newton-Raphson method with subsidiary function.
   *
   * \tparam SubsidiaryType Type of the subsidiary function.
   * \param subsidiaryFunction Subsidiary function to be solved.
      * \param x Where the solution should be stored.
   * \param subsidiaryArgs Additional arguments for the subsidiary function.
   * \param dxPredictor Predictor for the solution increment (default is NoPredictor).
   * \return Information about the solution process.
   */
  template <typename SubsidiaryType>
  [[nodiscard(
      "The solve method returns information of the solution process. You should store this information and check if "
      "it was successful")]] NonLinearSolverInformation
  solve(Domain& x,SubsidiaryType&& subsidiaryFunction, SubsidiaryArgs& subsidiaryArgs) {
    this->notify(NonLinearSolverMessages::INIT);

    /// Initializations
    Ikarus::NonLinearSolverInformation solverInformation;
    solverInformation.success = true;

    auto& lambda = x.parameter();

    /// Determine Fext0
    /// It is assumed that Fext = Fext0 * lambda such that dRdlambda = Fext0
    /// Generalization for Fext0 = Fext0(lambda) is not implemented

    auto lambdaDummy = lambda;
    auto DDummy      = x;
    x.setZero();
    lambda = 1.0;
    decltype(auto) rx = residualFunction_(x);
    const auto Fext0 = (-rx).eval(); // dRdlambda(lambda)
    lambda           = lambdaDummy;
    x                = DDummy;

    rx = residualFunction_(x);
    decltype(auto) Ax = jacobianFunction_(x);

    Eigen::VectorXd deltaD;
    deltaD.resizeLike(rx);
    deltaD.setZero();

    subsidiaryArgs.dfdDD.resizeLike(Fext0);

    subsidiaryFunction(subsidiaryArgs);
    auto rNorm = sqrt(rx.dot(rx));
    decltype(rNorm) dNorm;
    int iter{0};
    if constexpr (isLinearSolver)
      linearSolver_.analyzePattern(Ax);

    Eigen::MatrixX2<double> residual2d, sol2d;

    /// Iterative solving scheme
    while (rNorm > settings_.tol && iter < settings_.maxIter) {
      this->notify(NonLinearSolverMessages::ITERATION_STARTED);

      /// Two-step solving procedure
      residual2d.resize(rx.rows(), 2);
      residual2d << -rx, Fext0;
      sol2d.resize(rx.rows(), 2);

      if constexpr (isLinearSolver) {
        linearSolver_.factorize(Ax);
        linearSolver_.solve(sol2d, residual2d);
      } else {
        sol2d = linearSolver(residual2d, Ax);
      }

      subsidiaryFunction(subsidiaryArgs);

      const double deltalambda = (-subsidiaryArgs.f - subsidiaryArgs.dfdDD.dot(sol2d.col(0))) /
                                 (subsidiaryArgs.dfdDD.dot(sol2d.col(1)) + subsidiaryArgs.dfdDlambda);
      deltaD = sol2d.col(0) + deltalambda * sol2d.col(1);

      updateFunction_(x, deltaD);
      updateFunction_(subsidiaryArgs.DD, deltaD);

      lambda += deltalambda;
      subsidiaryArgs.Dlambda += deltalambda;

      dNorm = sqrt(deltaD.dot(deltaD) + deltalambda * deltalambda);
      rx = residualFunction_(x);
      Ax = jacobianFunction_(x);
      rNorm = sqrt(rx.dot(rx) + subsidiaryArgs.f * subsidiaryArgs.f);

      this->notify(NonLinearSolverMessages::SOLUTION_CHANGED, static_cast<double>(lambda));
      this->notify(NonLinearSolverMessages::CORRECTIONNORM_UPDATED, static_cast<double>(dNorm));
      this->notify(NonLinearSolverMessages::RESIDUALNORM_UPDATED, static_cast<double>(rNorm));
      this->notify(NonLinearSolverMessages::ITERATION_ENDED);

      ++iter;
    }

    if (iter == settings_.maxIter)
      solverInformation.success = false;
    solverInformation.iterations     = iter;
    solverInformation.residualNorm   = rNorm;
    solverInformation.correctionNorm = dNorm;
    if (solverInformation.success)
      this->notify(NonLinearSolverMessages::FINISHED_SUCESSFULLY, iter);

    return solverInformation;
  }

  /**
   * \brief Access the residual function.
   * \return Reference to the residual function.
   */
  auto& residual() { return residualFunction_; }

private:
  NonLinearOperator residualFunction_;
  typename NonLinearOperator::Derivative jacobianFunction_;
  CorrectionType correction_;
  LS linearSolver_;
  UF updateFunction_;
  Settings settings_;
};
/**
 * \brief Function to create a NewtonRaphson with subsidiary function solver instance.
 * \tparam NLO Type of the nonlinear operator to solve.
 * \tparam LS Type of the linear solver used internally (default is SolverDefault).
 * \tparam UF Type of the update function (default is UpdateDefault).
 * \param nonLinearOperator Nonlinear operator to solve.
 * \param linearSolver Linear solver used internally (default is SolverDefault).
 * \param updateFunction Update function (default is UpdateDefault).
 * \return Shared pointer to the NewtonRaphson solver instance.
 */
template <typename NLO, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
auto makeNewtonRaphsonWithSubsidiaryFunction(const NLO& nonLinearOperator, LS&& linearSolver = {},
                                             UF&& updateFunction = {}) {
  return std::make_shared<NewtonRaphsonWithSubsidiaryFunction<NLO, LS, UF>>(
      nonLinearOperator, std::forward<LS>(linearSolver), std::move(updateFunction));
}

template <typename NLO, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
NewtonRaphsonWithSubsidiaryFunction(const NLO& nonLinearOperator, LS&& linearSolver = {}, UF&& updateFunction = {})
    -> NewtonRaphsonWithSubsidiaryFunction<NLO, std::remove_cvref_t<LS>, std::remove_cvref_t<UF>>;

} // namespace Ikarus
