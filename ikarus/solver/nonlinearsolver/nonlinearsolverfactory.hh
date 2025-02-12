// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverfactory.hh
 * \brief Contains the generic NonlinearSolverFactory class.
 */

#pragma once

#include <utility>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/utils/defaultfunctions.hh>
#include <ikarus/utils/nonlinopfactory.hh>

namespace Ikarus {

/**
 * @brief A factory class for creating nonlinear solvers.
 *
 * This class is responsible for creating nonlinear solvers using the provided settings
 * and an assembler that satisfies the FlatAssembler concept.
 *
 * @tparam NLSSetting The type of the settings used for creating nonlinear solvers.
 */
template <typename NLSSetting>
struct NonlinearSolverFactory
{
  /**
   * @brief Constructs a NonlinearSolverFactory with the given settings.
   *
   * @param s The settings to be used by the factory.
   */
  NonlinearSolverFactory(NLSSetting s)
      : settings(s) {}

  NLSSetting settings;

  /**
   * @brief Creates a nonlinear solver using the provided assembler.
   *
   * The assembler must satisfy the FlatAssembler concept.
   *
   * @tparam Assembler The type of the assembler used for creating the nonlinear solver.
   * @param assembler The assembler to be used for creating the nonlinear solver.
   * @return The created nonlinear solver.
   *
   * @note The assembler's dBCOption is checked, and the appropriate update function
   *       is used based on whether the option is set to Reduced or not.
   */
  template <typename Assembler>
  requires Concepts::FlatAssembler<typename std::remove_cvref_t<Assembler>::element_type>
  auto create(Assembler&& assembler) const {
    auto nonLinOp         = NonLinearOperatorFactory::op(assembler);
    using SignatureT           = typename Dune::Functions::SignatureTraits<typename decltype(nonLinOp)::Signature>;
    using Domain           = typename SignatureT::Domain;
    using DerivativeSignatureT           = typename Dune::Functions::SignatureTraits<typename decltype(nonLinOp)::Derivative::Signature>;
    using CorrectionType           = typename DerivativeSignatureT::Range;
    std::function updateF = [assembler, setting = settings](Domain& a,
                                                            const CorrectionType& b) {
      if (assembler->dBCOption() == DBCOption::Reduced) {
        setting.updateFunction(a, assembler->createFullVector(b));
      } else
        setting.updateFunction(a, b);
    };
    auto settingsNew = settings.rebindUpdateFunction(std::move(updateF));
    return createNonlinearSolver(std::move(settingsNew), std::move(nonLinOp));
  }
};
}; // namespace Ikarus