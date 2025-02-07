// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controlroutinestate.hh
 * \brief State for all controlroutines
 */

#pragma once

#include <ikarus/utils/traits.hh>

namespace Ikarus {

// Forward declarations
#ifndef DOXYGEN
template <typename TypeListOne, typename TypeListTwo>
class NonLinearOperator;
#endif

/**
 * \brief State for control routines
 *
 * \tparam LoadParameter the type of the load parameter
 */
template <typename LoadParameter>
struct ControlRoutineState
{
  LoadParameter parameter;

  int loadStep{};
  double stepSize{};
};

namespace Impl {
  template <typename T>
  struct ControlRoutineStateFactory;

  template <typename NLO>
  requires traits::isSpecialization<NonLinearOperator, NLO>::value
  struct ControlRoutineStateFactory<NLO>
  {
  private:
    using LastParameter =
        const typename NLO::template ParameterValue<std::tuple_size_v<typename NLO::ParameterValues> - 1>&;

  public:
    using type = ControlRoutineState<LastParameter>;
  };
} // namespace Impl

/**
 * \brief Helper to deduce the correct types for ControlRoutineState
 *
 * \tparam NLO The nonlinear operator
 */
template <typename NLO>
using ControlRoutineStateType = Impl::ControlRoutineStateFactory<NLO>::type;

} // namespace Ikarus
