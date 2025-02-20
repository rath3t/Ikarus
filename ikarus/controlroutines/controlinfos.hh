// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file
 * \brief Defines the ControlInformation structure for storing control results.
 * \ingroup  controlroutines
 */

#pragma once
#include <vector>

#include <ikarus/solver/nonlinearsolver/solverinfos.hh>

namespace Ikarus {

/**
 * \struct ControlInformation
 * \brief Structure containing information about the control results.
 */
struct ControlInformation
{
  bool success{false}; ///< Flag indicating the success of the control.
  std::vector<Ikarus::NonLinearSolverInformation>
      solverInfos{};      ///< Vector containing information from nonlinear solvers.
  int totalIterations{0}; ///< Total number of iterations performed.
};

} // namespace Ikarus
