// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file trustregion.hh
 * \brief Implementation of the Trust-region method for solving nonlinear equations.
 */

#pragma once

#include <iosfwd>

#include <dune/common/float_cmp.hh>

#include <spdlog/spdlog.h>

#include <Eigen/Sparse>

#include <ikarus/linearalgebra/truncatedconjugategradient.hh>
#include <ikarus/solver/nonlinearsolver/helperfunctions.hh>
#include <ikarus/solver/nonlinearsolver/solverinfos.hh>
#include <ikarus/utils/defaultfunctions.hh>
#include <ikarus/utils/linearalgebrahelper.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>
#include <ikarus/utils/traits.hh>

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

struct TRSettings
{
  int verbosity    = 5;                                       ///< Verbosity level.
  double maxtime   = std::numeric_limits<double>::infinity(); ///< Maximum allowable time for solving.
  int minIter      = 3;                                       ///< Minimum number of iterations.
  int maxIter      = 1000;                                    ///< Maximum number of iterations.
  int debug        = 0;                                       ///< Debugging flag.
  double grad_tol  = 1e-6;                                    ///< Gradient tolerance.
  double corr_tol  = 1e-6;                                    ///< Correction tolerance.
  double rho_prime = 0.01;                                    ///< Rho prime value.
  bool useRand     = false;                                   ///< Flag for using random correction predictor.
  double rho_reg   = 1e6;                                     ///< Regularization value for rho.
  double Delta_bar = std::numeric_limits<double>::infinity(); ///< Maximum trust region radius.
  double Delta0    = 10;                                      ///< Initial trust region radius.
};

/**
 * \struct TrustRegionSettings
 * \brief Configuration settings for the TrustRegion solver.
 */

template <PreConditioner preConditioner = PreConditioner::IncompleteCholesky, typename UF = utils::UpdateDefault>
struct TrustRegionConfig
{
  static_assert(std::copy_constructible<UF>,
                " The update function should be copy constructable. If it is a lambda wrap it in a std::function");

  using UpdateFunction = UF;

  TRSettings parameters;
  static constexpr PreConditioner preConditionerType = preConditioner;
  UF updateFunction;
  template <typename UF2>
  auto rebindUpdateFunction(UF2&& updateFunction) const {
    TrustRegionConfig<preConditioner, UF2> settings{.parameters     = parameters,
                                                    .updateFunction = std::forward<UF2>(updateFunction)};
    return settings;
  }
};

template <typename NLO, PreConditioner preConditioner = PreConditioner::IncompleteCholesky,
          typename UF = utils::UpdateDefault>
class TrustRegion;

/**
 * \brief Function to create a trust region non-linear solver
 * \tparam NLO Type of the nonlinear operator to solve.
 * \tparam TRConfig Type of the nonlinear solver config.
 * \param config Config for the solver.
 * \param nonLinearOperator Nonlinear operator to solve.
 * \return Shared pointer to the TrustRegion solver instance.
 */
template <typename NLO, typename TRConfig>
requires traits::isSpecializationNonTypeAndTypes<TrustRegionConfig, std::remove_cvref_t<TRConfig>>::value
auto createNonlinearSolver(TRConfig&& config, NLO&& nonLinearOperator) {
  static constexpr PreConditioner preConditioner = std::remove_cvref_t<TRConfig>::preConditionerType;
  using UF                                       = std::remove_cvref_t<TRConfig>::UpdateFunction;
  static_assert(std::remove_cvref_t<NLO>::numberOfFunctions == 3,
                "The number of derivatives in the nonlinear operator have to be exactly 3.");
  auto solver = std::make_shared<TrustRegion<NLO, preConditioner, UF>>(nonLinearOperator,
                                                                       std::forward<TRConfig>(config).updateFunction);

  solver->setup(config.parameters);
  return solver;
}

/**
 * \enum StopReason
 * \brief Enumeration of reasons for stopping the TrustRegion solver.
 */
enum class StopReason
{
  gradientNormTolReached,
  correctionNormTolReached,
  maximumTimeReached,
  maximumIterationsReached,
  dontStop
};

/**
 * \struct AlgoInfo
 * \brief Additional information about the TrustRegion algorithm.
 */
struct AlgoInfo
{
  int consecutive_TRplus   = 0;                 ///< Consecutive trust region increases.
  int consecutive_TRminus  = 0;                 ///< Consecutive trust region decreases.
  int consecutive_Rejected = 0;                 ///< Consecutive rejected proposals.
  std::string stopReasonString;                 ///< String describing the stopping reason.
  std::string trstr;                            ///< Trust region change string (TR+, TR-).
  std::string accstr;                           ///< Acceptance string (REJ, acc).
  std::string randomPredictionString;           ///< Random prediction string.
  std::string cauchystr = "                  "; ///< Used Cauchy step string.
  bool acceptProposal;                          ///< Flag indicating whether the proposal is accepted.
  bool used_cauchy = false;                     ///< Flag indicating whether Cauchy point was used.
  StopReason stop{StopReason::dontStop};        ///< Stopping reason.
  std::string reasonString;                     ///< String describing the stopping reason.
};

/**
 * \struct Stats
 * \brief Information about the TrustRegion solver.
 */
struct Stats
{
  double gradNorm{1};
  double etaNorm{1};
  double time{0};
  double energy{0};
  double energyProposal{0};
  double rho{0};
  int outerIter{0};
  int innerIterSum{0};
};

/**
* \class TrustRegion
* \brief Trust Region solver for non-linear optimization problems.
* \details Refer to \cite trustregion for details of the algorithm.

* This code is heavily inspired by the trust-region implementation of
<a href="https://github.com/NicolasBoumal/manopt/blob/master/manopt/solvers/trustregions/trustregions.m">Manopt</a>.
* \ingroup solvers
* \tparam NLO Type of the nonlinear operator to solve.
* \tparam preConditioner Type of preconditioner to use (default is IncompleteCholesky).
* \tparam UF Type of the update function
*/
template <typename NLO, PreConditioner preConditioner, typename UF>
class TrustRegion : public IObservable<NonLinearSolverMessages>
{
public:
  using Settings  = TRSettings;                               ///< Type of the settings for the TrustRegion solver
  using ValueType = typename NLO::template ParameterValue<0>; ///< Type of the parameter vector of
                                                              ///< the nonlinear operator
  using CorrectionType = typename NLO::DerivativeType;        ///< Type of the correction of x += deltaX.
  using UpdateFunction = UF;                                  ///< Type of the update function.

  using NonLinearOperator = NLO; ///< Type of the non-linear operator

  using ScalarType = std::remove_cvref_t<typename NLO::template FunctionReturnType<0>>; ///< Type of the scalar
                                                                                        ///< cost

  using MatrixType = std::remove_cvref_t<typename NLO::template FunctionReturnType<2>>; ///< Type of the Hessian

  /**
   * \brief Constructs a TrustRegion solver instance.
   * \param nonLinearOperator Nonlinear operator to solve.
   * \param updateFunction Update function
   */
  template <typename UF2 = UF>
  explicit TrustRegion(const NLO& nonLinearOperator, UF2&& updateFunction = {})
      : nonLinearOperator_{nonLinearOperator},
        updateFunction_{std::forward<UF2>(updateFunction)},
        xOld_{this->nonLinearOperator().firstParameter()} {
    eta_.setZero(gradient().size());
    Heta_.setZero(gradient().size());
    truncatedConjugateGradient_.analyzePattern(hessian());
  }

  /**
   * \brief Sets up the TrustRegion solver with the provided settings and checks feasibility.
   * \param p_settings TrustRegionSettings containing the solver configuration.
   */
  void setup(const Settings& settings) {
    settings_ = settings;
    assert(settings_.rho_prime < 0.25 && "options.rho_prime must be strictly smaller than 1/4.");
    assert(settings_.Delta_bar > 0 && "options.Delta_bar must be positive.");
    assert(settings_.Delta0 > 0 && settings_.Delta0 < settings_.Delta_bar &&
           "options.Delta0 must be positive and smaller than Delta_bar.");
  }

#ifndef DOXYGEN
  struct NoPredictor
  {
  };
#endif
  /**
   * \brief Solves the nonlinear optimization problem using the TrustRegion algorithm.
   * \tparam SolutionType Type of the solution predictor (default is NoPredictor).
   * \param dxPredictor Solution predictor.
   * \return NonLinearSolverInformation containing information about the solver result.
   */
  template <typename SolutionType = NoPredictor>
  requires std::is_same_v<SolutionType, NoPredictor> || std::is_convertible_v<SolutionType, CorrectionType>
  NonLinearSolverInformation solve(const SolutionType& dxPredictor = NoPredictor{}) {
    this->notify(NonLinearSolverMessages::INIT);
    stats_ = Stats{};
    info_  = AlgoInfo{};

    NonLinearSolverInformation solverInformation;
    nonLinearOperator().updateAll();
    stats_.energy   = energy();
    auto& x         = nonLinearOperator().firstParameter();
    xOld_           = x;
    stats_.gradNorm = norm(gradient());
    if constexpr (not std::is_same_v<SolutionType, NoPredictor>)
      updateFunction(x, dxPredictor);
    truncatedConjugateGradient_.analyzePattern(hessian());

    innerInfo_.Delta = settings_.Delta0;
    spdlog::info(
        "        | iter | inner_i |   rho |   energy | energy_p | energy_inc |  norm(g) |    Delta | norm(corr) | "
        "InnerBreakReason");
    spdlog::info("{:-^143}", "-");
    while (not stoppingCriterion()) {
      this->notify(NonLinearSolverMessages::ITERATION_STARTED);
      if (settings_.useRand) {
        if (stats_.outerIter == 0) {
          eta_.setRandom();
          while (eta_.dot(hessian() * eta_) > innerInfo_.Delta * innerInfo_.Delta)
            eta_ *= eps_; // eps is sqrt(sqrt(maschine-precision))
        } else
          eta_.setZero();
      } else
        eta_.setZero();

      solveInnerProblem();
      stats_.innerIterSum += innerInfo_.numInnerIter;

      info_.stopReasonString = tcg_stop_reason_[static_cast<int>(innerInfo_.stop_tCG)];
      Heta_                  = hessian() * eta_;
      if (settings_.useRand and stats_.outerIter == 0) {
        info_.used_cauchy            = false;
        info_.randomPredictionString = " Used Random correction predictor";
        info_.cauchystr              = "                  ";
        double tauC;
        // Check the curvature
        const Eigen::VectorXd Hg = hessian() * gradient();
        const auto g_Hg          = (gradient().dot(Hg));
        if (g_Hg <= 0)
          tauC = 1;
        else
          tauC = std::min(Dune::power(stats_.gradNorm, 3) / (innerInfo_.Delta * g_Hg), 1.0);

        // generate the Cauchy point.
        const Eigen::VectorXd etaC  = -tauC * innerInfo_.Delta / stats_.gradNorm * gradient();
        const Eigen::VectorXd HetaC = -tauC * innerInfo_.Delta / stats_.gradNorm * Hg;

        const double mdle  = stats_.energy + gradient().dot(eta_) + .5 * Heta_.dot(eta_);
        const double mdlec = stats_.energy + gradient().dot(etaC) + .5 * HetaC.dot(etaC);
        if (mdlec < mdle && stats_.outerIter == 0) {
          eta_              = etaC;
          Heta_             = HetaC;
          info_.used_cauchy = true;

          info_.cauchystr = " USED CAUCHY POINT!";
        }
      } else
        info_.cauchystr = "                  ";

      stats_.etaNorm = eta_.norm();

      updateFunction_(x, eta_);

      // Calculate energy of our proposed update step
      nonLinearOperator().template update<0>();
      stats_.energyProposal = energy();

      // Will we accept the proposal or not?
      // Check the performance of the quadratic model against the actual energy.
      auto rhonum = stats_.energy - stats_.energyProposal;
      auto rhoden = -eta_.dot(gradient() + 0.5 * Heta_);

      /*  Close to convergence the proposed energy and the real energy almost coincide.
       *  Therefore, the performance check of our model becomes ill-conditioned
       *  The regularisation fixes this */
      const auto rhoReg = std::max(1.0, abs(stats_.energy)) * eps_ * settings_.rho_reg;
      rhonum            = rhonum + rhoReg;
      rhoden            = rhoden + rhoReg;

      const bool modelDecreased = rhoden > 0.0;

      if (!modelDecreased)
        info_.stopReasonString.append(", model did not decrease");

      stats_.rho = rhonum / rhoden;
      stats_.rho = stats_.rho < 0.0 ? -1.0 : stats_.rho; // move rho to the domain [-1.0,inf]

      info_.trstr = "   ";

      // measure if energy decreased
      const bool energyDecreased = Dune::FloatCmp::ge(stats_.energy - stats_.energyProposal, -1e-12);

      // If the model behaves badly or if the energy increased we reduce the trust region size
      if (stats_.rho < 1e-4 || not modelDecreased || std::isnan(stats_.rho) || not energyDecreased) {
        info_.trstr = "TR-";
        innerInfo_.Delta /= 4.0;
        info_.consecutive_TRplus = 0;
        info_.consecutive_TRminus++;
        if (info_.consecutive_TRminus >= 5 && settings_.verbosity >= 1) {
          info_.consecutive_TRminus = -std::numeric_limits<int>::infinity();
          spdlog::info(" +++ Detected many consecutive TR- (radius decreases).");
          spdlog::info(" +++ Consider decreasing options.Delta_bar by an order of magnitude.");
        }

      } else if (stats_.rho > 0.99 && (innerInfo_.stop_tCG == Eigen::TCGStopReason::negativeCurvature ||
                                       innerInfo_.stop_tCG == Eigen::TCGStopReason::exceededTrustRegion)) {
        info_.trstr               = "TR+";
        innerInfo_.Delta          = std::min(3.5 * innerInfo_.Delta, settings_.Delta_bar);
        info_.consecutive_TRminus = 0;
        info_.consecutive_TRplus++;
        if (info_.consecutive_TRplus >= 5 && settings_.verbosity >= 1) {
          info_.consecutive_TRplus = -std::numeric_limits<int>::infinity();
          spdlog::info(" +++ Detected many consecutive TR+ (radius increases)");
          spdlog::info(" +++ Consider increasing options.Delta_bar by an order of magnitude");

        } else {
          info_.consecutive_TRplus  = 0;
          info_.consecutive_TRminus = 0;
        }
      }

      if (modelDecreased && stats_.rho > settings_.rho_prime && energyDecreased) {
        if (stats_.energyProposal > stats_.energy)
          spdlog::info(
              "Energy function increased by {} (step size: {}). Since this is small we accept the step and hope for "
              "convergence of the gradient norm.",
              stats_.energyProposal - stats_.energy, stats_.etaNorm);

        info_.acceptProposal       = true;
        info_.accstr               = "acc";
        info_.consecutive_Rejected = 0;
      } else {
        info_.acceptProposal = false;
        info_.accstr         = "REJ";

        if (info_.consecutive_Rejected >= 5)
          innerInfo_.Delta /= 2;
        else
          innerInfo_.Delta = std::min(innerInfo_.Delta, stats_.etaNorm / 2.0);
        ++info_.consecutive_Rejected;
      }

      stats_.outerIter++;

      if (settings_.verbosity == 1)
        logState();

      info_.randomPredictionString = "";

      updateStates(nonLinearOperator().assembler(), eta_);
      if (info_.acceptProposal) {
        stats_.energy = stats_.energyProposal;
        nonLinearOperator_.updateAll();
        xOld_ = x;
        this->notify(NonLinearSolverMessages::CORRECTIONNORM_UPDATED, stats_.etaNorm);
        this->notify(NonLinearSolverMessages::RESIDUALNORM_UPDATED, stats_.gradNorm);
        this->notify(NonLinearSolverMessages::SOLUTION_CHANGED);
      } else {
        x = xOld_;
        eta_.setZero();
      }
      nonLinearOperator_.updateAll();
      stats_.gradNorm = gradient().norm();
      this->notify(NonLinearSolverMessages::ITERATION_ENDED);
    }
    spdlog::info("{}", info_.reasonString);
    spdlog::info("Total iterations: {} Total CG Iterations: {}", stats_.outerIter, stats_.innerIterSum);

    solverInformation.success =
        (info_.stop == StopReason::correctionNormTolReached) or (info_.stop == StopReason::gradientNormTolReached);

    solverInformation.iterations   = stats_.outerIter;
    solverInformation.residualNorm = stats_.gradNorm;
    if (solverInformation.success)
      this->notify(NonLinearSolverMessages::FINISHED_SUCESSFULLY, solverInformation.iterations);
    return solverInformation;
  }
  /**
   * \brief Access the nonlinear operator.
   * \return Reference to the nonlinear operator.
   */
  auto& nonLinearOperator() { return nonLinearOperator_; }

private:
  void logState() const {
    spdlog::info(
        "{:>3s} {:>3s} {:>6d} {:>9d}  {:>6.2f}  {:>9.2e}  {:>9.2e}  {:>11.2e}  {:>9.2e}  {:>9.2e}  {:>11.2e}   "
        "{:<73}",
        info_.accstr, info_.trstr, stats_.outerIter, innerInfo_.numInnerIter, stats_.rho, stats_.energy,
        stats_.energyProposal, stats_.energyProposal - stats_.energy, stats_.gradNorm, innerInfo_.Delta, stats_.etaNorm,
        info_.stopReasonString + info_.cauchystr + info_.randomPredictionString);
  }

  void logFinalState() {
    spdlog::info("{:>3s} {:>3s} {:>6d} {:>9d}  {: ^6}  {: ^9}  {: ^9}  {: ^11}  {:>9.2e}  {: ^9}  {: ^11}   {:<73}",
                 info_.accstr, info_.trstr, stats_.outerIter, innerInfo_.numInnerIter, " ", " ", " ", " ",
                 stats_.gradNorm, " ", " ", info_.stopReasonString + info_.cauchystr + info_.randomPredictionString);
  }

  inline const auto& energy() { return nonLinearOperator().value(); }
  inline const auto& gradient() { return nonLinearOperator().derivative(); }
  inline const auto& hessian() { return nonLinearOperator().secondDerivative(); }

  bool stoppingCriterion() {
    std::ostringstream stream;
    /** Gradient correction tolerance reached  */
    if (stats_.gradNorm < settings_.grad_tol && stats_.outerIter != 0) {
      logFinalState();
      spdlog::info("CONVERGENCE:  Energy: {:1.16e}    norm(gradient): {:1.16e}", nonLinearOperator().value(),
                   stats_.gradNorm);
      stream << "Gradient norm tolerance reached; options.tolerance = " << settings_.grad_tol;

      info_.reasonString = stream.str();

      info_.stop = StopReason::gradientNormTolReached;
      return true;
    } else if (stats_.etaNorm < settings_.corr_tol && stats_.outerIter != 0) {
      logFinalState();
      spdlog::info("CONVERGENCE:  Energy: {:1.16e}    norm(correction): {:1.16e}", nonLinearOperator().value(),
                   stats_.etaNorm);
      stream << "Displacement norm tolerance reached;  = " << settings_.corr_tol << "." << std::endl;

      info_.reasonString = stream.str();
      info_.stop         = StopReason::correctionNormTolReached;
      return true;
    }

    /** Maximum Time reached  */
    if (stats_.time >= settings_.maxtime) {
      logFinalState();
      stream << "Max time exceeded; options.maxtime = " << settings_.maxtime << ".";
      info_.reasonString = stream.str();
      info_.stop         = StopReason::maximumTimeReached;
      return true;
    }

    /** Maximum Iterations reached  */
    if (stats_.outerIter >= settings_.maxIter) {
      logFinalState();
      stream << "Max iteration count reached; options.maxiter = " << settings_.maxIter << ".";
      info_.reasonString = stream.str();
      info_.stop         = StopReason::maximumIterationsReached;
      return true;
    }
    return false;
  }

  void solveInnerProblem() {
    truncatedConjugateGradient_.setInfo(innerInfo_);
    int attempts = 0;
    truncatedConjugateGradient_.factorize(hessian());
    // If the preconditioner is IncompleteCholesky the factorization may fail if we have negative diagonal entries and
    // the initial shift is too small. Therefore, if the factorization fails we increase the initial shift by a factor
    // of 5.
    if constexpr (preConditioner == PreConditioner::IncompleteCholesky) {
      while (truncatedConjugateGradient_.info() != Eigen::Success) {
        choleskyInitialShift_ *= 5;
        truncatedConjugateGradient_.preconditioner().setInitialShift(choleskyInitialShift_);
        truncatedConjugateGradient_.factorize(hessian());
        if (attempts > 5)
          DUNE_THROW(Dune::MathError, "Factorization of preconditioner failed!");
        ++attempts;
      }
      if (truncatedConjugateGradient_.info() == Eigen::Success)
        choleskyInitialShift_ = 1e-3;
    }
    eta_       = truncatedConjugateGradient_.solveWithGuess(-gradient(), eta_);
    innerInfo_ = truncatedConjugateGradient_.getInfo();
  }

  NLO nonLinearOperator_;
  UpdateFunction updateFunction_;
  typename NLO::template ParameterValue<0> xOld_;
  CorrectionType eta_;
  CorrectionType Heta_;
  Settings settings_;
  AlgoInfo info_;
  double choleskyInitialShift_ = 1e-3;
  Eigen::TCGInfo<double> innerInfo_;
  Stats stats_;
  static constexpr double eps_ = 0.0001220703125; // 0.0001220703125 is sqrt(sqrt(maschine-precision))
  std::array<std::string, 6> tcg_stop_reason_{
      {"negative curvature", "exceeded trust region", "reached target residual-kappa (linear)",
       "reached target residual-theta (superlinear)", "maximum inner iterations", "model increased"}
  };

  using PreConditionerType =
      std::conditional_t<preConditioner == PreConditioner::IdentityPreconditioner, Eigen::IdentityPreconditioner,
                         std::conditional_t<preConditioner == PreConditioner::DiagonalPreconditioner,
                                            typename Eigen::DiagonalPreconditioner<ScalarType>,
                                            typename Eigen::IncompleteCholesky<ScalarType>>>;
  Eigen::TruncatedConjugateGradient<MatrixType, Eigen::Lower | Eigen::Upper, PreConditionerType>
      truncatedConjugateGradient_;
};

/**
 * \brief Creates an instance of the TrustRegion solver.
 *
 * \tparam NLO Type of the nonlinear operator to solve.
 * \tparam preConditioner Type of the preconditioner used internally (default is IncompleteCholesky).
 * \tparam UF Type of the update function (default is UpdateDefault).
 * \param nonLinearOperator Nonlinear operator to solve.
 * \param updateFunction Update function (default is UpdateDefault).
 * \return Shared pointer to the TrustRegion solver instance.
 */
template <typename NLO, PreConditioner preConditioner = PreConditioner::IncompleteCholesky,
          typename UF = utils::UpdateDefault>
auto makeTrustRegion(const NLO& nonLinearOperator, UF&& updateFunction = {}) {
  return std::make_shared<TrustRegion<NLO, preConditioner, UF>>(nonLinearOperator, updateFunction);
}

template <typename NLO, PreConditioner preConditioner = PreConditioner::IncompleteCholesky,
          typename UF2 = utils::UpdateDefault>
TrustRegion(const NLO& nonLinearOperator,
            UF2&& updateFunction = {}) -> TrustRegion<NLO, preConditioner, std::remove_cvref_t<UF2>>;

} // namespace Ikarus
