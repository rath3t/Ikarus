// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file derivativetraits.hh
 * \brief Contains derivative traits for common vbalue and derivative relations used for
 * makeDifferentiableFunctionFromCallables
 */

#pragma once

#include <dune/functions/common/differentiablefunctionfromcallables.hh>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <ikarus/utils/traits.hh>

namespace Ikarus {

namespace Impl {
  template <typename... Args>
  struct Functions;
} // namespace Impl

enum class FEParameter;
enum class FESolutions;
template <FESolutions sol, FEParameter para, typename SV, typename PM>
class FERequirements;

/**
 * \brief Represents a NonLinearOperator class for handling nonlinear operators.
 * \ingroup utils
 * \tparam TypeListOne The type list for the first set of functions.
 * \tparam TypeListTwo The type list for the second set of functions.
 */
template <typename TypeListOne, typename TypeListTwo>
class DerivativeTraitsFromCallables
{
public:
  DerivativeTraitsFromCallables([[maybe_unused]] const TypeListOne& derivativesFunctions,
                                [[maybe_unused]] const TypeListTwo& args) {
    static_assert(!sizeof(TypeListOne),
                  "This type should not be instantiated. check that your arguments satisfies the template below");
  }
};

template <typename... DerivativeArgs, typename  Arg>
struct DerivativeTraitsFromCallables<Impl::Functions<DerivativeArgs...>, Arg>
{
  /**
   * \brief Constructor for DerivativeTraitsFromCallables.
   *
   * \param derivativesFunctions The Functions object for derivative arguments.
   * \param parameterI The Parameter object for parameter arguments.
   */
   template<typename Arg2=Arg>
  DerivativeTraitsFromCallables(const Impl::Functions<DerivativeArgs...>& derivativesFunctions,
                                 Arg2&& parameterI) {}
  using Ranges = std::tuple<std::invoke_result_t<DerivativeArgs, Arg>...,Dune::Functions::InvalidRange>;
  using Domain = std::remove_cvref_t<Arg>;
  using Signatures = std::tuple<std::invoke_result_t<DerivativeArgs, Arg>(Arg)...,Dune::Functions::InvalidRange(Arg)>;
  using RawSignatures = std::tuple<std::remove_cvref_t<std::invoke_result_t<DerivativeArgs, Arg>>(std::remove_cvref_t<Arg>)...,Dune::Functions::InvalidRange(Arg)>;

  template <int I>
  using Signature = std::tuple_element_t<I/*(I <std::tuple_size_v<Signatures> ? I : std::tuple_size_v<Signatures>-1)*/, Signatures>;

  template <int I>
  using Range = std::tuple_element_t<I, Ranges>;

  template <typename Signature>
  struct DerivativeTraits
  {
  private:
    static constexpr int indexOfSignatureImpl = traits::Index<Signature, RawSignatures>::value+1;
    static constexpr int indexOfSignature = indexOfSignatureImpl;//<std::tuple_size_v<Ranges>? indexOfSignatureImpl: std::tuple_size_v<Ranges>-1;

  public:
    using Range = std::tuple_element_t<indexOfSignature, Ranges>;
  };
};
template <typename... DerivativeArgs, typename  Arg>
  DerivativeTraitsFromCallables(const Impl::Functions<DerivativeArgs...>& derivativesFunctions,
                                 Arg&& parameterI) -> DerivativeTraitsFromCallables<Impl::Functions<DerivativeArgs...>, Arg>;
} // namespace Ikarus
