// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file derivativetraits.hh
 * \brief Contains derivative traits for common vbalue and derivative relations used for makeDifferentiableFunctionFromCallables
 */

#pragma once


#include <dune/functions/common/differentiablefunctionfromcallables.hh>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Ikarus {

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

template <typename... DerivativeArgs, typename... ParameterArgs>
struct DerivativeTraitsFromCallables<Impl::Functions<DerivativeArgs...>, Impl::Parameter<ParameterArgs...>>
{
    /**
   * \brief Constructor for DerivativeTraitsFromCallables.
   *
   * \param derivativesFunctions The Functions object for derivative arguments.
   * \param parameterI The Parameter object for parameter arguments.
   */
  template <typename U = void>
  requires(not std::is_rvalue_reference_v<DerivativeArgs> and ...)
  explicit DerivativeTraitsFromCallables(const Impl::Functions<DerivativeArgs...>& derivativesFunctions,
                             const Impl::Parameter<ParameterArgs...>& parameterI)
      {}

    using FunctionReturnValuesWrapper = std::tuple<std::invoke_result_t<DerivativeArgs, ParameterArgs...>...>;


  using Signatures = std::tuple<std::invoke_result_t<DerivativeArgs, ParameterArgs...>(ParameterArgs...)...>;
  template<int I>
  using Signature = std::tuple_element_t<I,Signatures>;

template<typename Signature>
    struct DerivativeTraits
    {
        static constexpr int indexOfSignature = traits::Index<Signature,Signatures>::value;
        using Range = std::conditional_t<indexOfSignature<std::tuple_size_v<FunctionReturnValuesWrapper> , std::tuple_element_t<indexOfSignature,FunctionReturnValuesWrapper>,Dune::Functions::InvalidRange>;
    };
  
  

};
        template<class Signature>
    struct DerivativeTraitsDense
    {
      typedef Dune::Functions::InvalidRange Range;
    };

template<typename K> requires std::floating_point<K>
struct DerivativeTraitsDense< K(K) >
{
  typedef double Range;
};

template<typename K, int n>
struct DerivativeTraitsDense<K(Eigen::Vector<K,n>)>
{
  typedef Eigen::Vector<K,n> Range;
};

template <typename K, FESolutions sol, FEParameter para, typename SV, typename PM>
struct DerivativeTraitsDense<K(FERequirements<sol,para,SV,PM>)>
{
  typedef FERequirements<sol,para,SV,PM>::SolutionVectorType Range;
};

template<typename K,  int n,  int m>
struct DerivativeTraitsDense<Eigen::Vector<K,m>(Eigen::Vector<K,n> )>
{
  typedef Eigen::Matrix<K,m,n> Range;
};

template <typename K,  int n,FESolutions sol, FEParameter para, typename SV, typename PM>
struct DerivativeTraitsDense<Eigen::Vector<K,n>(FERequirements<sol,para,SV,PM>)>
{
  typedef Eigen::Matrix<K,n,n> Range;
};

struct DerivativeTraitsDenseD
{
  template<typename D>
  using Traits= DerivativeTraitsDense<D>;
};


      template<class Signature>
    struct DerivativeTraitsSparse: DerivativeTraitsDense<Signature>
    {
    };

template<typename K,  int n,  int m>
struct DerivativeTraitsSparse<Eigen::Vector<K,m>(Eigen::Vector<K,n>)>
{
  typedef Eigen::SparseMatrix<K> Range;
};

template <typename K,  int n,FESolutions sol, FEParameter para, typename SV, typename PM>
struct DerivativeTraitsSparse<Eigen::Vector<K,n>(FERequirements<sol,para,SV,PM>)>
{
  typedef Eigen::SparseMatrix<K> Range;
};

struct DerivativeTraitsSparseD
{
  template<typename D>
  using Traits= DerivativeTraitsSparse<D>;
};


} // namespace Ikarus
