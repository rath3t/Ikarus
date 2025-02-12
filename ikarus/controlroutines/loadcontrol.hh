// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file loadcontrol.hh
 * \brief Defines the LoadControl class.
 */

#pragma once

#include <memory>

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
  LoadControl(const std::shared_ptr<NLS>& nonLinearSolver, int loadSteps, const std::array<double, 2>& tbeginEnd)
      : nonLinearSolver_{nonLinearSolver},
        loadSteps_{loadSteps},
        parameterBegin_{tbeginEnd[0]},
        parameterEnd_{tbeginEnd[1]},
        stepSize_{(parameterEnd_ - parameterBegin_) / loadSteps_} {  }

  /**
   * \brief Executes the LoadControl routine.
   *  \tparam Domain Type of the domain.
    * \param x The solution.
   * \return ControlInformation structure containing information about the control results.
   */
   template<typename Domain>
  ControlInformation run(Domain& x);

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
