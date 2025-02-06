// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controlroutinefactory.hh
 * \brief Factory for controlroutines
 */

#pragma once
#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/concepts.hh>

namespace Ikarus {

/**
 * \brief A factory class for creating control routines
 *
 * \details This class is responsible for creating control routines using the provided config and an assembler.
 * The assembler is used to register the finite elements to the broadcaster messages from the solver and the
 * control routine.
 *
 */

struct ControlRoutineFactory
{
  /**
   * \brief Creates the control routine and registering the elements to the broadcaster messages from the solver and
   * the control routine.
   *
   * \tparam NLS the type of the nonlinear solver used by the control routine
   * \tparam Assembler the type of the assembler. It has to be a pointer and satisfy the FlatAssembler conncept
   * \tparam CRConfig The Type of the config
   * \param config the config
   * \param nonlinearSolver the nonlinear solver
   * \param assembler the assembler
   * \return The control routine
   */
  template <typename CRConfig, typename NLS, typename Assembler>
  requires Concepts::FlatAssembler<typename std::remove_cvref_t<Assembler>::element_type>
  static auto create(const CRConfig& config, NLS&& nonlinearSolver, Assembler&& assembler) {
    auto cr = createControlRoutine(std::move(config), std::forward<NLS>(nonlinearSolver));

    for (auto& fe : assembler->finiteElements()) {
      fe.template subscribeTo<NonLinearSolverMessages>(nonlinearSolver);
      fe.template subscribeTo<ControlMessages>(cr);
    }
    return cr;
  }
};

} // namespace Ikarus
