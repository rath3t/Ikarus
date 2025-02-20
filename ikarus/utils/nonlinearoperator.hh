// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearoperator.hh
 * \brief Provides a NonLinearOperator class for handling nonlinear operators.
 *
 */

#pragma once
#include <tuple>
#include <type_traits>

#include "dune/functions/common/differentiablefunctionfromcallables.hh"
#include <dune/common/hybridutilities.hh>

#include <ikarus/utils/derivativetraits.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {

namespace Impl {

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
} // namespace Impl


/**
 * \brief Creates a Functions object.
 *
 * \tparam Args The argument types.
 * \param args The tuple of arguments.
 * \return auto The Functions object.
 */
template <typename... Args>
auto functions(Args&&... args) {
  return Impl::Functions<std::remove_cvref_t<Args>...>{std::forward_as_tuple(std::forward<Args>(args)...)};
}
#ifndef DOXYGEN
template <class Signature, template <class> class DerivativeTraits, class... F>
class NonLinearOperator;


#endif

template <class Range, class D, template <class> class DerivativeTraits, class F>
class NonLinearOperator<Range(D), DerivativeTraits, F>
    : private Dune::Functions::DifferentiableFunctionFromCallables<Range(D), DerivativeTraits, F>
{

  using Base = Dune::Functions::DifferentiableFunctionFromCallables<Range(D), DerivativeTraits, F>;
  public:
  using Domain = D;
    static constexpr auto nDerivatives = 0;

  template<class FFF, Dune::disableCopyMove<NonLinearOperator, FFF> = 0>
  NonLinearOperator(FFF&& f) :
    Base(std::forward<FFF>(f))
  {}
  using Traits = DerivativeTraitsFromCallables<Impl::Functions<F>, Domain>;

  using Derivative = NonLinearOperator<typename Traits::template Signature<1>, DerivativeTraits>;

  Range operator() (const D& x) const
  {
    return Base::operator()(x );
  }

    /**
   * \brief Get derivative of NonLinearOperator
   *
   */
   friend Derivative derivative(const NonLinearOperator& t)
   {
     return derivative(static_cast<const Base&>(t));
   }
};

/**
 * \brief NonLinearOperator is a class taking linear algebra function and their arguments.
 * The function are assumed to be derivatives of each other w.r.t. the first parameter
 *
 * \tparam DerivativeArgs The types of derivative arguments.
 * \tparam ParameterArgs The types of parameter arguments.
 */
template <class Range, class D, template <class> class DerivativeTraits, class F,class... FF>
class NonLinearOperator<Range(D), DerivativeTraits, F,FF...>
    : private Dune::Functions::DifferentiableFunctionFromCallables<Range(D), DerivativeTraits, F,FF...>
{

  using Base = Dune::Functions::DifferentiableFunctionFromCallables<Range(D), DerivativeTraits, F,FF...>;
  public:
    using Domain = std::remove_cvref_t<D>;
    static constexpr auto nDerivatives = sizeof...(FF);

  template<class... FFF>
  NonLinearOperator(FFF&&... f) :
    Base(std::forward<FFF>(f)...)
  {}
  using Traits = DerivativeTraitsFromCallables<Impl::Functions<F,FF...>, D>;


  using Derivative = NonLinearOperator<typename Traits::template Signature<1>, DerivativeTraits, FF...>;

  Range operator() (const D& x) const
  {
    return Base::operator()(x );
  }

    /**
   * \brief Get derivative of NonLinearOperator
   *
   */
   friend Derivative derivative(const NonLinearOperator& t)
   {
     auto df= derivative(static_cast<const Base&>(t));
      return Derivative(df);
   }
};

template <typename... DerivativeArgs, typename Arg>
auto makeNonLinearOperator(const Impl::Functions<DerivativeArgs...>& derivativesFunctions,
                            Arg&& parameter) {
  DerivativeTraitsFromCallables t(derivativesFunctions, std::forward<Arg>(parameter));
  using DerivTraits = decltype(t);
  auto la           = []<typename... F>(F&&... f) {
    return Ikarus::NonLinearOperator<typename DerivTraits::template Signature<0>,
                                               DerivTraits::template DerivativeTraits, std::remove_cvref_t<F>...>(std::forward<F>(f)...);
  };
  return std::apply(la, derivativesFunctions.args);
}

} // namespace Ikarus
