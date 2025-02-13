// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearoperator.hh
 * \brief Provides a NonLinearOperator class for handling nonlinear operators.
 *
 */

#pragma once
#include <functional>

#include <dune/common/hybridutilities.hh>
#include <dune/functions/common/differentiablefunctionfromcallables.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {

namespace Impl {
  /**
   * \brief Helper function to apply a function and remove reference wrappers.
   *
   * \tparam F The function type.
   * \tparam Tuple The tuple type.
   * \tparam I Index sequence for tuple elements.
   * \param f The function to be applied.
   * \param t The tuple of arguments.
   * \return constexpr decltype(auto) The result after applying the function and removing reference wrappers.
   */
  template <class F, class Tuple, std::size_t... I>
  constexpr decltype(auto) applyAndRemoveRefererenceWrapper(F&& f, Tuple&& t, std::index_sequence<I...>) {
    return std::invoke(std::forward<F>(f),
                       std::get<I>(std::forward<Tuple>(t)).get()...); //.get gets the impl type of std::referenceWrapper
  }

  /**
   * \brief Helper function to apply a function and remove reference wrappers.
   *
   * \tparam F The function type.
   * \tparam Tuple The tuple type.
   * \param f The function to be applied.
   * \param t The tuple of arguments.
   * \return constexpr decltype(auto) The result after applying the function and removing reference wrappers.
   */
  template <class F, class Tuple>
  constexpr decltype(auto) applyAndRemoveReferenceWrapper(F&& f, Tuple&& t) {
    return applyAndRemoveRefererenceWrapper(
        std::forward<F>(f), std::forward<Tuple>(t),
        std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<Tuple>>>{});
  }


  /**
   * \brief Helper function to make a tuple of values and reference wrappers.
   *
   * \tparam Pars The tuple type.
   * \tparam Tuple The tuple of arguments.
   * \tparam I Index sequence for tuple elements.
   * \param t The tuple of arguments.
   * \param p The tuple of parameters.
   * \return constexpr decltype(auto) The resulting tuple of values and reference wrappers.
   */
  template <class Pars, class Tuple, std::size_t... I>
  constexpr decltype(auto) makeTupleOfValuesAndReferences(Tuple&& t, Pars&& p, std::index_sequence<I...>) {
    return std::make_tuple(
        forwardasReferenceWrapperIfIsReference(applyAndRemoveReferenceWrapper(std::get<I>(t), p))...);
  }

  /**
   * \brief Represents a tuple of functions.
   *
   * \tparam Args The argument types.
   */
  template <typename... Args>
  struct Functions
  {
    std::tuple<Args...> args;
  };

  /**
   * \brief Represents a tuple of parameters.
   *
   * \tparam Args The argument types.
   */
  template <typename... Args>
  struct Parameter
  {
    std::tuple<Args...> args;
  };

} // namespace Impl

/**
 * \brief Creates a Parameter object.
 *
 * \tparam Args The argument types.
 * \param args The tuple of arguments.
 * \return auto The Parameter object.
 */
template <typename... Args>
auto parameter(Args&&... args) {
  return Impl::Parameter<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
}

/**
 * \brief Creates a Functions object.
 *
 * \tparam Args The argument types.
 * \param args The tuple of arguments.
 * \return auto The Functions object.
 */
template <typename... Args>
auto functions(Args&&... args) {
  return Impl::Functions<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
}

/**
 * \brief Initializes the results for functions and parameters.
 *
 * \tparam DerivativeArgs The types of derivative arguments.
 * \tparam ParameterArgs The types of parameter arguments.
 * \param derivativesFunctions The Functions object for derivative arguments.
 * \param parameter The Parameter object for parameter arguments.
 * \return auto The initialized results.
 */
template <typename... DerivativeArgs, typename... ParameterArgs>
auto initResults(const std::tuple<DerivativeArgs...>& derivativesFunctions,
                 const std::tuple<ParameterArgs...>& parameter) {
  return Impl::makeTupleOfValuesAndReferences(
      derivativesFunctions, parameter,
      std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<decltype(derivativesFunctions)>>>{});
}

/**
 * \brief NonLinearOperator is a class taking linear algebra function and their arguments.
 * The function are assumed to be derivatives of each other w.r.t. the first parameter
 *
 * \tparam DerivativeArgs The types of derivative arguments.
 * \tparam ParameterArgs The types of parameter arguments.
 */
 template<class Range, class Domain, template<class> class DerivativeTraits, class F>
class NonLinearOperator : public Dune::Functions::DifferentiableFunctionFromCallables<Range(Domain), DerivativeTraits, F>
{
using Base= Dune::Functions::DifferentiableFunctionFromCallables<Range(Domain), DerivativeTraits, F>;
using Base::Base;



};

} // namespace Ikarus
