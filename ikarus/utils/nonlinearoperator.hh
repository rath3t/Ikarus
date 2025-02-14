// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearoperator.hh
 * \brief Provides a NonLinearOperator class for handling nonlinear operators.
 *
 */

#pragma once
#include <tuple>

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
#ifndef DOXYGEN
template <class Signature, template <class> class DerivativeTraits, class... F>
class NonLinearOperator;


#endif

template <class Range, class Domain, template <class> class DerivativeTraits, class F>
class NonLinearOperator<Range(Domain), DerivativeTraits, F>
    : private Dune::Functions::DifferentiableFunctionFromCallables<Range(Domain), DerivativeTraits, F>
{
  
  using Base = Dune::Functions::DifferentiableFunctionFromCallables<Range(Domain), DerivativeTraits, F>;
  
  public:
  template<class FFF, Dune::disableCopyMove<NonLinearOperator, FFF> = 0>
  NonLinearOperator(FFF&& f) :
    Base(std::forward<FFF>(f))
  {}
  using Traits = DerivativeTraitsFromCallables<Impl::Functions<F>, Impl::Parameter<Domain>>;

  using Derivative = NonLinearOperator<typename Traits::template Signature<1>, DerivativeTraits>;

  Range operator() (const Domain& x) const
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
template <class Range, class Domain, template <class> class DerivativeTraits, class F,class... FF>
class NonLinearOperator<Range(Domain), DerivativeTraits, F,FF...>
    : private Dune::Functions::DifferentiableFunctionFromCallables<Range(Domain), DerivativeTraits, F,FF...>
{
  
  using Base = Dune::Functions::DifferentiableFunctionFromCallables<Range(Domain), DerivativeTraits, F,FF...>;
  
  public:
  template<class... FFF>
  NonLinearOperator(FFF&&... f) :
    Base(std::forward<FFF>(f)...)
  {}
  using Traits = DerivativeTraitsFromCallables<Impl::Functions<F,FF...>, Impl::Parameter<Domain>>;

  using Derivative = NonLinearOperator<typename Traits::template Signature<1>, DerivativeTraits, FF...>;

  Range operator() (const Domain& x) const
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

template <typename... DerivativeArgs, typename... ParameterArgs>
auto makeNonLinearOperator(const Impl::Functions<DerivativeArgs...>& derivativesFunctions,
                           const Impl::Parameter<ParameterArgs...>& parameters) {
  DerivativeTraitsFromCallables t(derivativesFunctions, parameters);
  using DerivTraits = decltype(t);
  auto la           = []<typename... F>(F&&... f) {
    return Ikarus::NonLinearOperator<typename DerivTraits::template Signature<0>,
                                               DerivTraits::template DerivativeTraits, F...>(std::forward<F>(f)...);
  };
  return std::apply(la, derivativesFunctions.args);
}

} // namespace Ikarus
