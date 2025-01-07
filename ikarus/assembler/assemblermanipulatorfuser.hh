// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file assemblermanipulatorfuser.hh
 * \brief Defines a decorator for the assemblers that helps to manipulate the assembled quantities
 */

#pragma once
#include <ikarus/assembler/assemblermanipulatorbuildingblocks.hh>

namespace Ikarus {

namespace Impl {
  template <template <typename, typename, typename, typename> typename AInterface,
            template <typename Wrapper, typename Assembler> typename AImplementation>
  struct AssemblerInterfaceHelper
  {
    template <typename WA, typename A, typename ReturnValueType>
    using Interface = AInterface<WA, typename A::FEContainer, typename A::DirichletValuesType, ReturnValueType>;
    template <typename WA, typename A>
    using Implementation = AImplementation<WA, A>;
  };

} // namespace Impl
/**
 * \class AssemblerManipulator
 * \brief The AssemblerManipulator defines a decorator for the assemblers that
 * helps to manipulate the assembled quantities.
 * \ingroup assembler
 *
 * \tparam A Type of the assembler.
 */

template <Concepts::FlatAssembler A, typename... Ass>
class AssemblerManipulator;

#define BASECLASSMEMBERFUNCTION(func, A)         \
  template <typename... Args>                    \
  decltype(auto) base_##func(Args... args) {     \
    return A::func(std::forward<Args>(args)...); \
  }

/**
 * \class AssemblerManipulator
 * \brief Specialization for handling scalar assembly manipulations using callback functions
 * \ingroup assembler
 *
 * \tparam A Type of the assembler.
 * \tparam ScalarAss Scalar assembler type.
 */
template <Concepts::FlatAssembler A, typename ScalarAss>
class AssemblerManipulator<A, ScalarAss>
    : public ScalarAss::template Interface<AssemblerManipulator<A, ScalarAss>, A, typename A::ScalarType>,
      public ScalarAss::template Implementation<AssemblerManipulator<A, ScalarAss>, A>,
      private A
{
public:
  using WrappedAssembler = A;

protected:
  using ScalarAssemblerImpl      = ScalarAss::template Implementation<AssemblerManipulator, WrappedAssembler>;
  using ScalarAssemblerInterface = ScalarAss::template Interface<AssemblerManipulator, A, typename A::ScalarType>;

  friend ScalarAssemblerImpl;
  friend ScalarAssemblerInterface;

public:
  // typdefs from FlatAssemblerBase
  using typename WrappedAssembler::AffordanceCollectionType;
  using typename WrappedAssembler::Basis;
  using typename WrappedAssembler::DirichletValuesType;
  using typename WrappedAssembler::FEContainer;
  using typename WrappedAssembler::FERequirement;
  using typename WrappedAssembler::ScalarType;
  using typename WrappedAssembler::SizeType;

  // functions from FlatAssemblerBase
  using WrappedAssembler::bind;
  using WrappedAssembler::bound;
  using WrappedAssembler::boundToAffordanceCollection;
  using WrappedAssembler::boundToDBCOption;
  using WrappedAssembler::boundToRequirement;
  using WrappedAssembler::createFullVector;

  using WrappedAssembler::affordanceCollection;
  using WrappedAssembler::constraintsBelow;
  using WrappedAssembler::dBCOption;
  using WrappedAssembler::estimateOfConnectivity;
  using WrappedAssembler::gridView;
  using WrappedAssembler::isConstrained;
  using WrappedAssembler::reducedSize;
  using WrappedAssembler::requirement;
  using WrappedAssembler::size;

  // functions from implementations
  using ScalarAssemblerImpl::bind;

private:
  using ScalarAssemblerImpl::getScalarImpl;

public:
  // functions from interfaces of the new manipulator assembler
  using ScalarAssemblerInterface::scalar;

  using CallBackTypes =
      std::tuple<typename ScalarAss::template Implementation<AssemblerManipulator<A, ScalarAss>, A>::FunctionType>;

  /**
   * \brief Constructor that forwards arguments to the base assembler
   *
   * \tparam Args Variadic template parameter for constructor arguments
   * \param args Arguments to be forwarded to the base assembler
   */
  template <typename... Args>
  requires(not std::is_same_v<std::remove_cvref_t<std::tuple_element_t<0, std::tuple<Args...>>>, AssemblerManipulator>)
  explicit AssemblerManipulator(Args&&... args)
      : WrappedAssembler(std::forward<Args>(args)...) {}

private:
  BASECLASSMEMBERFUNCTION(getScalarImpl, A)

  A& base() { return *this; }
  const A& base() const { return *this; }
};

/**
 * \class AssemblerManipulator
 * \brief Specialization for handling scalar and vector assembly manipulations using callback functions
 * \ingroup assembler
 *
 * \tparam A Type of the assembler.
 * \tparam ScalarAss Scalar assembler type.
 * \tparam VectorAss Vector assembler type.
 */
template <Concepts::FlatAssembler A, typename ScalarAss, typename VectorAss>
class AssemblerManipulator<A, ScalarAss, VectorAss>
    : public ScalarAss::template Interface<AssemblerManipulator<A, ScalarAss, VectorAss>, A, typename A::ScalarType>,
      public ScalarAss::template Implementation<AssemblerManipulator<A, ScalarAss, VectorAss>, A>,
      public VectorAss::template Interface<AssemblerManipulator<A, ScalarAss, VectorAss>, A, typename A::VectorType>,
      public VectorAss::template Implementation<AssemblerManipulator<A, ScalarAss, VectorAss>, A>,
      private A
{
public:
  using WrappedAssembler = A;

protected:
  using ScalarAssemblerImpl      = ScalarAss::template Implementation<AssemblerManipulator, WrappedAssembler>;
  using VectorAssemblerImpl      = VectorAss::template Implementation<AssemblerManipulator, WrappedAssembler>;
  using ScalarAssemblerInterface = ScalarAss::template Interface<AssemblerManipulator, A, typename A::ScalarType>;
  using VectorAssemblerInterface = VectorAss::template Interface<AssemblerManipulator, A, typename A::VectorType>;

  friend ScalarAssemblerImpl;
  friend VectorAssemblerImpl;
  friend ScalarAssemblerInterface;
  friend VectorAssemblerInterface;

public:
  // typdefs from FlatAssemblerBase
  using typename WrappedAssembler::AffordanceCollectionType;
  using typename WrappedAssembler::Basis;
  using typename WrappedAssembler::DirichletValuesType;
  using typename WrappedAssembler::FEContainer;
  using typename WrappedAssembler::FERequirement;
  using typename WrappedAssembler::ScalarType;
  using typename WrappedAssembler::SizeType;
  using typename WrappedAssembler::VectorType;

  // functions from FlatAssemblerBase
  using WrappedAssembler::bind;
  using WrappedAssembler::bound;
  using WrappedAssembler::boundToAffordanceCollection;
  using WrappedAssembler::boundToDBCOption;
  using WrappedAssembler::boundToRequirement;
  using WrappedAssembler::createFullVector;

  using WrappedAssembler::affordanceCollection;
  using WrappedAssembler::constraintsBelow;
  using WrappedAssembler::dBCOption;
  using WrappedAssembler::estimateOfConnectivity;
  using WrappedAssembler::gridView;
  using WrappedAssembler::isConstrained;
  using WrappedAssembler::reducedSize;
  using WrappedAssembler::requirement;
  using WrappedAssembler::size;

  // functions from implementations
  using ScalarAssemblerImpl::bind;
  using VectorAssemblerImpl::bind;

private:
  using ScalarAssemblerImpl::getScalarImpl;
  using VectorAssemblerImpl::getRawVectorImpl;
  using VectorAssemblerImpl::getReducedVectorImpl;
  using VectorAssemblerImpl::getVectorImpl;

public:
  // functions from interfaces of the new manipulator assembler
  using ScalarAssemblerInterface::scalar;
  using VectorAssemblerInterface::vector;

  using CallBackTypes =
      std::tuple<typename ScalarAssemblerImpl::FunctionType, typename VectorAssemblerImpl::FunctionType>;

  /**
   * \brief Constructor that forwards arguments to the base assembler
   *
   * \tparam Args Variadic template parameter for constructor arguments
   * \param args Arguments to be forwarded to the base assembler
   */
  template <typename... Args>
  requires(not std::is_same_v<std::remove_cvref_t<std::tuple_element_t<0, std::tuple<Args...>>>, AssemblerManipulator>)
  explicit AssemblerManipulator(Args&&... args)
      : WrappedAssembler(std::forward<Args>(args)...) {}

private:
  BASECLASSMEMBERFUNCTION(getScalarImpl, WrappedAssembler)
  BASECLASSMEMBERFUNCTION(getRawVectorImpl, WrappedAssembler)
  BASECLASSMEMBERFUNCTION(getReducedVectorImpl, WrappedAssembler)
  BASECLASSMEMBERFUNCTION(getVectorImpl, WrappedAssembler)

  WrappedAssembler& base() { return *this; }
  const WrappedAssembler& base() const { return *this; }
};

/**
 * \class AssemblerManipulator
 * \brief Specialization for handling scalar and vector assembly manipulations using callback functions
 * \ingroup assembler
 *
 * \tparam A Type of the assembler.
 * \tparam ScalarAss Scalar assembler type.
 * \tparam VectorAss Vector assembler type.
 * \tparam MatrixAss Matrix assembler type.
 */
template <Concepts::FlatAssembler A, typename ScalarAss, typename VectorAss, typename MatrixAss>
class AssemblerManipulator<A, ScalarAss, VectorAss, MatrixAss>
    : public ScalarAss::template Interface<AssemblerManipulator<A, ScalarAss, VectorAss, MatrixAss>, A,
                                           typename A::ScalarType>,
      public ScalarAss::template Implementation<AssemblerManipulator<A, ScalarAss, VectorAss, MatrixAss>, A>,
      public VectorAss::template Interface<AssemblerManipulator<A, ScalarAss, VectorAss, MatrixAss>, A,
                                           typename A::VectorType>,
      public VectorAss::template Implementation<AssemblerManipulator<A, ScalarAss, VectorAss, MatrixAss>, A>,
      public MatrixAss::template Interface<AssemblerManipulator<A, ScalarAss, VectorAss, MatrixAss>, A,
                                           typename A::MatrixType>,
      public MatrixAss::template Implementation<AssemblerManipulator<A, ScalarAss, VectorAss, MatrixAss>, A>,
      private A
{
public:
  using WrappedAssembler = A;

protected:
  using ScalarAssemblerImpl      = ScalarAss::template Implementation<AssemblerManipulator, A>;
  using VectorAssemblerImpl      = VectorAss::template Implementation<AssemblerManipulator, A>;
  using MatrixAssemblerImpl      = MatrixAss::template Implementation<AssemblerManipulator, A>;
  using ScalarAssemblerInterface = ScalarAss::template Interface<AssemblerManipulator, A, typename A::ScalarType>;
  using VectorAssemblerInterface = VectorAss::template Interface<AssemblerManipulator, A, typename A::VectorType>;
  using MatrixAssemblerInterface = MatrixAss::template Interface<AssemblerManipulator, A, typename A::MatrixType>;

  friend ScalarAssemblerImpl;
  friend VectorAssemblerImpl;
  friend MatrixAssemblerImpl;
  friend MatrixAssemblerInterface;
  friend ScalarAssemblerInterface;
  friend VectorAssemblerInterface;

public:
  // typdefs from FlatAssemblerBase
  using typename WrappedAssembler::AffordanceCollectionType;
  using typename WrappedAssembler::Basis;
  using typename WrappedAssembler::DirichletValuesType;
  using typename WrappedAssembler::FEContainer;
  using typename WrappedAssembler::FERequirement;
  using typename WrappedAssembler::MatrixType;
  using typename WrappedAssembler::ScalarType;
  using typename WrappedAssembler::SizeType;
  using typename WrappedAssembler::VectorType;

  // functions from FlatAssemblerBase
  using WrappedAssembler::bind;
  using WrappedAssembler::bound;
  using WrappedAssembler::boundToAffordanceCollection;
  using WrappedAssembler::boundToDBCOption;
  using WrappedAssembler::boundToRequirement;
  using WrappedAssembler::createFullVector;

  using WrappedAssembler::affordanceCollection;
  using WrappedAssembler::constraintsBelow;
  using WrappedAssembler::dBCOption;
  using WrappedAssembler::estimateOfConnectivity;
  using WrappedAssembler::gridView;
  using WrappedAssembler::isConstrained;
  using WrappedAssembler::reducedSize;
  using WrappedAssembler::requirement;
  using WrappedAssembler::size;

  // functions from implementations
  using MatrixAssemblerImpl::bind;
  using ScalarAssemblerImpl::bind;
  using VectorAssemblerImpl::bind;

private:
  using MatrixAssemblerImpl::getMatrixImpl;
  using MatrixAssemblerImpl::getRawMatrixImpl;
  using MatrixAssemblerImpl::getReducedMatrixImpl;
  using ScalarAssemblerImpl::getScalarImpl;
  using VectorAssemblerImpl::getRawVectorImpl;
  using VectorAssemblerImpl::getReducedVectorImpl;
  using VectorAssemblerImpl::getVectorImpl;

public:
  // functions from interfaces of the new manipulator assembler
  using MatrixAssemblerInterface::matrix;
  using ScalarAssemblerInterface::scalar;
  using VectorAssemblerInterface::vector;

  using CallBackTypes =
      std::tuple<typename ScalarAssemblerImpl::FunctionType, typename VectorAssemblerImpl::FunctionType,
                 typename MatrixAssemblerImpl::FunctionType>;

  /**
   * \brief Constructor that forwards arguments to the base assembler
   *
   * \tparam Args Variadic template parameter for constructor arguments
   * \param args Arguments to be forwarded to the base assembler
   */
  template <typename... Args>
  requires(not std::is_same_v<std::remove_cvref_t<std::tuple_element_t<0, std::tuple<Args...>>>, AssemblerManipulator>)
  explicit AssemblerManipulator(Args&&... args)
      : WrappedAssembler(std::forward<Args>(args)...) {}

private:
  BASECLASSMEMBERFUNCTION(getScalarImpl, A)
  BASECLASSMEMBERFUNCTION(getRawVectorImpl, A)
  BASECLASSMEMBERFUNCTION(getReducedVectorImpl, A)
  BASECLASSMEMBERFUNCTION(getVectorImpl, A)
  BASECLASSMEMBERFUNCTION(getRawMatrixImpl, A)

  BASECLASSMEMBERFUNCTION(getReducedMatrixImpl, A)

  BASECLASSMEMBERFUNCTION(getMatrixImpl, A)

  A& base() { return *this; }
  const A& base() const { return *this; }
};

/**
 * \brief Creates an AssemblerManipulator instance based on the type of the provided assembler.
 *
 * This function constructs an appropriate AssemblerManipulator for the given assembler object,
 * utilizing different interface and implementation templates based on the assembler's capabilities.
 *
 * \tparam A Type of the assembler.
 * \param a An assembler object to be wrapped by the AssemblerManipulator.
 * \return An instance of AssemblerManipulator wrapping the provided assembler.
 */
template <typename A>
auto makeAssemblerManipulator(A&& a) {
  constexpr bool scal = Concepts::ScalarFlatAssembler<std::remove_cvref_t<A>>;
  constexpr bool vec  = Concepts::VectorFlatAssembler<std::remove_cvref_t<A>>;
  constexpr bool mat  = Concepts::MatrixFlatAssembler<std::remove_cvref_t<A>>;
  if constexpr (scal && vec && mat)
    return std::make_shared<
        AssemblerManipulator<std::remove_cvref_t<A>, Impl::AssemblerInterfaceHelper<ScalarAssembler, ScalarManipulator>,
                             Impl::AssemblerInterfaceHelper<VectorAssembler, VectorManipulator>,
                             Impl::AssemblerInterfaceHelper<MatrixAssembler, MatrixManipulator>>>(std::forward<A>(a));
  else if constexpr (scal && vec)
    return std::make_shared<
        AssemblerManipulator<std::remove_cvref_t<A>, Impl::AssemblerInterfaceHelper<ScalarAssembler, ScalarManipulator>,
                             Impl::AssemblerInterfaceHelper<VectorAssembler, VectorManipulator>>>(std::forward<A>(a));
  else if constexpr (scal)
    return std::make_shared<AssemblerManipulator<std::remove_cvref_t<A>,
                                                 Impl::AssemblerInterfaceHelper<ScalarAssembler, ScalarManipulator>>>(
        std::forward<A>(a));
}
} // namespace Ikarus
