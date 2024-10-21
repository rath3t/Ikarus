// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file helperfunctions.hh
 * \brief Implementation of the helper functions used by the nonlinear solvers.
 */

#pragma once

#include <ikarus/utils/concepts.hh>
namespace Ikarus {
/**
 * \brief A helper function that calls the updateState() function of the finite elements to update the internal
 * variables.
 * \tparam Assembler Type of the assembler.
 * \param assembler Shared pointer to the underlying assembler
 * \param correction A correction vector (for example, the displacement increment) used to update the internal state
 * variables of the finite elements.
 */
template <typename Assembler>
void updateStates(std::shared_ptr<Assembler>& assembler, const auto& correction) {
  if constexpr (Concepts::FlatAssembler<Assembler>) {
    auto fes = assembler->finiteElements();
    auto req = assembler->requirement();
    for (auto& fe : fes)
      updateState(fe, req, correction);
  } else {
    /* Dummy branch used if std::is_same_v<Assembler, Impl::NoAssembler> */
  }
}
} // namespace Ikarus