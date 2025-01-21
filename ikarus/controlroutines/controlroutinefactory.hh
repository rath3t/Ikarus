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
 * The assembler is used to register the finite elements to the broadcaster messenges from the solver and the
 * control routine.
 *
 * \tparam CRConfig The Type of the config
 */
template <typename CRConfig>
struct ControlRoutineFactory
{
  /**
   * \brief Construct a new Control Routine Factory with the given config
   *
   * \param config the config used by the factory
   */
  ControlRoutineFactory(const CRConfig& config)
      : config_(config) {}

  /**
   * \brief Creates the control routine and registering the elements to the broadcaster messenges from the solver and
   * the control routine.
   *
   * \tparam NLS the type of the nonlinear solver used by the control routine
   * \tparam Assembler the type of the assembler. It has to be a pointer and statisfy the FlatAssembler conncept
   * \param nonlinearSolver the nonlinenar sovler
   * \param assembler the assembler
   * \return The control routine
   */
  template <typename NLS, typename Assembler>
  requires Concepts::FlatAssembler<typename std::remove_cvref_t<Assembler>::element_type>
  auto create(NLS&& nonlinearSolver, Assembler&& assembler) const {
    auto cr = createControlRoutine(std::move(config_), std::forward<NLS>(nonlinearSolver));

    for (auto& fe : assembler->finiteElements()) {
      fe.template subscribeTo<NonLinearSolverMessages>(nonlinearSolver);
      fe.template subscribeTo<ControlMessages>(cr);
    }

    return cr;
  }

private:
  CRConfig config_;
};

} // namespace Ikarus
