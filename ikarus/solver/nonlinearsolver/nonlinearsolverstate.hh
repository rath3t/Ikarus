// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverstate.hh
 * \brief State for all nonlinear solvers
 */

#pragma once

#include <ikarus/utils/traits.hh>

namespace Ikarus {

// Forward declarations
template <typename TypeListOne, typename TypeListTwo>
class NonLinearOperator;

enum class FESolutions;
enum class FEParameter;
template <FESolutions sol, FEParameter para, typename SV, typename PM>
class FERequirements;

template <typename P1, typename SolType = P1>
struct NonlinearSolverState
{
  P1 correction;
  SolType solution;

  double rNorm{};
  double dNorm{};
  int iteration{};
};

template <typename T>
struct NonlinearSolverStateFactory;

template <typename NLO>
requires traits::isSpecialization<NonLinearOperator, NLO>::value
struct NonlinearSolverStateFactory<NLO>
{
  using type = NonlinearSolverState<const typename NLO::ValueType&, const typename NLO::template ParameterValue<0>&>;
};

template <typename FER>
requires traits::isSpecializationNonTypeNonTypeAndTypes<FERequirements, FER>::value
struct NonlinearSolverStateFactory<FER>
{
  using type = NonlinearSolverState<const typename FER::SolutionVectorType&>;
};

template <typename T>
using NonlinearSolverStateType = NonlinearSolverStateFactory<T>::type;

} // namespace Ikarus
