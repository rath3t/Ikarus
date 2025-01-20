// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file loadcontrol.hh
 * \brief Defines the LoadControl class.
 */

#pragma once

#include <memory>

#include <dune/common/hybridutilities.hh>

#include <ikarus/controlroutines/controlinfos.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

/**
 * \class LoadControl
 * \brief The LoadControl control routine increases the last parameter of a nonlinear operator and calls a nonlinear
 * solver.
 *   \ingroup controlroutines
 * This class represents the LoadControl control routine. It increments the last parameter of a nonlinear operator
 * and utilizes a nonlinear solver, such as Newton's method, to solve the resulting system at each step.
 *
 * \tparam NLS Type of the nonlinear solver used in the control routine.
 */
template <typename NLS>
class LoadControl : public IObservable<ControlMessages>
{
public:
  /** \brief The name of the LoadControl method. */
  constexpr auto name() const { return std::string("Load Control Method"); }

  /**
   * \brief Constructor for LoadControl.
   *
   * \param nonLinearSolver Shared pointer to the nonlinear solver.
   * \param loadSteps Number of load steps in the control routine.
   * \param tbeginEnd Array representing the range of load parameters [tbegin, tend].
   */
  template <typename Assembler>
  LoadControl(const std::shared_ptr<NLS>& nonLinearSolver, int loadSteps, const std::array<double, 2>& tbeginEnd,
              std::shared_ptr<Assembler>& assembler)
      : nonLinearSolver_{nonLinearSolver},
        loadSteps_{loadSteps},
        parameterBegin_{tbeginEnd[0]},
        parameterEnd_{tbeginEnd[1]},
        stepSize_{(parameterEnd_ - parameterBegin_) / loadSteps_} {
    static_assert(
        requires {
          nonLinearSolver_->nonLinearOperator().lastParameter() = 0.0;
          nonLinearSolver_->nonLinearOperator().lastParameter() += 0.0;
        }, "The last parameter (load factor) must be assignable and incrementable with a double!");

    // register FEs to listen to NR messages
    if constexpr (not std::is_same_v<Assembler, Impl::NoAssembler>) {
      for (auto& fe : assembler->finiteElements()) {
        Dune::Hybrid::forEach(fe.getSubsciptions(),
                              [&](auto& subscription) { fe.subscribe(*nonLinearSolver_, subscription); });
      };
    }
  }

  /**
   * \brief Executes the LoadControl routine.
   *
   * \return ControlInformation structure containing information about the control results.
   */
  ControlInformation run();

  /* \brief returns the nonlinear solver */
  NLS& nonlinearSolver() { return *nonLinearSolver_; }

private:
  std::shared_ptr<NLS> nonLinearSolver_;
  int loadSteps_;
  double parameterBegin_;
  double parameterEnd_;
  double stepSize_;
};

} // namespace Ikarus

#include <ikarus/controlroutines/loadcontrol.inl>
