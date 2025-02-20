// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file concepts.hh
 * \brief Several concepts
 */

#pragma once

#include <concepts>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <autodiff/forward/dual/dual.hpp>

#include "ikarus/assembler/dirichletbcenforcement.hh"
#include "ikarus/finiteelements/mechanics/materials/tags.hh"
#include <ikarus/utils/traits.hh>

namespace Eigen {
template <typename Derived>
struct EigenBase;
}

namespace Ikarus {

template <typename Derived>
auto transpose(const Eigen::EigenBase<Derived>& A);
namespace Concepts {

  /**
   * \concept FlatInterLeavedBasis
   * \brief Concept to check if a basis uses FlatInterleaved indexing strategy.
   *
   * This concept checks if the given Basis type uses FlatInterleaved indexing strategy.
   *
   * \tparam Basis The basis type.
   */
  template <typename Basis>
  concept FlatInterLeavedBasis = requires {
    std::is_same_v<typename Basis::PreBasis::IndexMergingStrategy, Dune::Functions::BasisFactory::FlatInterleaved>;
  };

  namespace Impl {
    template <template <typename, int, typename> class U, typename T>
    struct LagrangeNodeHelper : std::false_type
    {
    };
    template <template <typename, int, typename> class U, typename GV, int k, typename R>
    struct LagrangeNodeHelper<U, U<GV, k, R>> : std::true_type
    {
    };

    template <template <typename, int, typename> class U, typename T, int k>
    struct LagrangeNodeHelperOfOrder : std::false_type
    {
    };
    template <template <typename, int, typename> class U, typename GV, int k, typename R>
    struct LagrangeNodeHelperOfOrder<U, U<GV, k, R>, k> : std::true_type
    {
    };
  } // namespace Impl

  /**
   * \concept LagrangeNode
   * \brief Concept to check if a node in a basis tree is a Lagrangian node.
   *
   *
   * \tparam N The node.
   */
  template <typename N>
  concept LagrangeNode = Impl::LagrangeNodeHelper<Dune::Functions::LagrangeNode, N>::value;

  /**
   * \concept LagrangeNode
   * \brief Concept to check if a node in a basis tree is a Lagrangian node with specific order.
   *
   *
   * \tparam N The node.
   */
  template <typename N, int order>
  concept LagrangeNodeOfOrder = Impl::LagrangeNodeHelperOfOrder<Dune::Functions::LagrangeNode, N, order>::value;

  /**
   * \concept FlatLexicographicBasis
   * \brief Concept to check if a basis uses FlatLexicographic indexing strategy.
   *
   * This concept checks if the given Basis type uses FlatLexicographic indexing strategy.
   *
   * \tparam B The basis type.
   */
  template <typename B>
  concept FlatLexicographicBasis = requires {
    std::is_same_v<typename B::PreBasis::IndexMergingStrategy, Dune::Functions::BasisFactory::FlatLexicographic>;
  };

  /**
   * \concept FlatIndexBasis
   * \brief Concept to check if a basis uses FlatIndex indexing strategy.
   *
   * This concept checks if the given Basis type uses either FlatLexicographic or FlatInterleaved indexing strategy.
   *
   * \tparam B The basis type.
   */
  template <typename B>
  concept FlatIndexBasis = FlatLexicographicBasis<B> or FlatInterLeavedBasis<B>;

  /**
   * \concept BlockedInterLeavedBasis
   * \brief Concept to check if a basis uses BlockedInterleaved indexing strategy.
   *
   * This concept checks if the given Basis type uses BlockedInterleaved indexing strategy.
   *
   * \tparam Basis The basis type.
   */
  template <typename Basis>
  concept BlockedInterLeavedBasis = requires {
    std::is_same_v<typename Basis::PreBasis::IndexMergingStrategy, Dune::Functions::BasisFactory::BlockedInterleaved>;
  };

  /**
   * \concept BlockedLexicographicBasis
   * \brief Concept to check if a basis uses BlockedLexicographic indexing strategy.
   *
   * This concept checks if the given Basis type uses BlockedLexicographic indexing strategy.
   *
   * \tparam B The basis type.
   */
  template <typename B>
  concept BlockedLexicographicBasis = requires {
    std::is_same_v<typename B::PreBasis::IndexMergingStrategy, Dune::Functions::BasisFactory::BlockedLexicographic>;
  };

  /**
   * \concept DuneLocalBasis
   * \brief Concept to check if a local basis is a duneLocalBasis.
   * \tparam DLB The dune local basis type .
   */
  template <typename DLB>
  concept DuneLocalBasis = requires(DLB& duneLocalBasis) {
    typename DLB::Traits::RangeType;
    typename DLB::Traits::JacobianType;
    DLB::Traits::dimDomain;
    typename DLB::Traits::DomainType;

    typename DLB::Traits::DomainFieldType;
    typename DLB::Traits::RangeFieldType;

    duneLocalBasis.evaluateFunction(std::declval<typename DLB::Traits::DomainType>(),
                                    std::declval<std::vector<typename DLB::Traits::RangeType>&>());
    duneLocalBasis.evaluateJacobian(std::declval<typename DLB::Traits::DomainType>(),
                                    std::declval<std::vector<typename DLB::Traits::JacobianType>&>());
  };

  /**
   * \concept BlockedIndexBasis
   * \brief Concept to check if a basis uses either BlockedLexicographic or BlockedInterleaved indexing strategy.
   *
   * This concept checks if the given Basis type uses either BlockedLexicographic or BlockedInterleaved indexing
   * strategy.
   *
   * \tparam B The basis type.
   */
  template <typename B>
  concept BlockedIndexBasis = BlockedLexicographicBasis<B> or BlockedInterLeavedBasis<B>;

  /**
   * \concept PathFollowingStrategy
   * \brief Concept defining the requirements for a path-following strategy.
   * \tparam PF Type representing the path-following strategy.
   * \tparam NLO Type representing the non-linear operator.
   * \tparam SA Type representing the subsidiary arguments.
   */
  template <typename PF, typename NLO, typename SA>
  concept PathFollowingStrategy = requires(PF pft, NLO nop, SA args) {
    { pft(args) } -> std::same_as<void>;
    { pft.initialPrediction(nop, args) } -> std::same_as<void>;
    { pft.intermediatePrediction(nop, args) } -> std::same_as<void>;
  };

  /**
   * \concept AdaptiveStepSizingStrategy
   * \brief Concept to check if a type implements all the needed functions to be an adaptive step sizing method.
   *
   * \tparam ASS The adaptive step sizing type.
   * \tparam NLSI The non-linear solver information type.
   * \tparam SA The subsidiary arguments type.
   */
  template <typename ASS, typename NLSI, typename SA, typename NonLinearOperator>
  concept AdaptiveStepSizingStrategy = requires(ASS adaptiveStepSizing, NLSI info, SA args, NonLinearOperator nop) {
    { adaptiveStepSizing(info, args, nop) } -> std::same_as<void>;
    { adaptiveStepSizing.targetIterations() } -> std::same_as<int>;
    { adaptiveStepSizing.setTargetIterations(std::declval<int>()) } -> std::same_as<void>;
  };

  /**
   * \concept LinearSolverCheck
   * \brief Concept to check if a linear solver implements all the needed functions for given vector and matrix types.
   *
   * \tparam LS The linear solver type.
   * \tparam M The matrix type.
   * \tparam V The vector type.
   */
  template <typename LS, typename M, typename V>
  concept LinearSolverCheck = requires(LS& linearSolver, M& A, V& vec) {
    linearSolver.analyzePattern(A);
    linearSolver.factorize(A);
    linearSolver.solve(vec, vec);
  };

  /**
   * \concept NonLinearSolverCheckForPathFollowing
   * \brief Concept to check if a non-linear solver with its non-linear operator satisfies  requirements for path
   * following.
   *
   * \tparam NLS The non-linear solver type.
   */
  template <typename NLS>
  concept NonLinearSolverCheckForPathFollowing = requires {
    std::tuple_size<typename NLS::NonLinearOperator::ParameterValues>::value == 2;
    not(std::is_same_v<typename NLS::NonLinearOperator::ValueType, double> and
        ((traits::isSpecializationTypeAndNonTypes<Eigen::Matrix,
                                                  typename NLS::NonLinearOperator::DerivativeType>::value) or
         (traits::isSpecializationTypeNonTypeAndType<Eigen::SparseMatrix,
                                                     typename NLS::NonLinearOperator::DerivativeType>::value)));
  };

  /**
   * \concept MultiplyAble
   * \brief Concept defining the requirements for types that support multiplication.
   * \tparam L Type of the left operand.
   * \tparam R Type of the right operand.
   * \details The concept specifies that instances of types `L` and `R` can be multiplied using the `*` operator.
   */
  template <typename L, typename R>
  concept MultiplyAble = requires(L x, R y) { x* y; };

  /**
   * \concept AddAble
   * \brief Concept defining the requirements for types that support addition.
   * \tparam L Type of the left operand.
   * \tparam R Type of the right operand.
   * \details The concept specifies that instances of types `L` and `R` can be added using the `+` operator.
   */
  template <typename L, typename R>
  concept AddAble = requires(L x, R y) { x + y; };

  /**
   * \concept SubstractAble
   * \brief Concept defining the requirements for types that support subtraction.
   * \tparam L Type of the left operand.
   * \tparam R Type of the right operand.
   * \details The concept specifies that instances of types `L` and `R` can be subtracted using the `-` operator.
   */
  template <typename L, typename R>
  concept SubstractAble = requires(L x, R y) { x - y; };

  /**
   * \concept MultiplyAssignAble
   * \brief Concept defining the requirements for types that support in-place multiplication.
   * \tparam L Type of the left operand.
   * \tparam R Type of the right operand.
   * \details The concept specifies that instances of type `L` can be multiplied in-place by instances of type `R`
   * using the `*=` operator.
   */
  template <typename L, typename R>
  concept MultiplyAssignAble = requires(L x, R y) { x *= y; };

  /**
   * \concept DivideAssignAble
   * \brief Concept defining the requirements for types that support in-place division.
   * \tparam L Type of the left operand.
   * \tparam R Type of the right operand.
   * \details The concept specifies that instances of type `L` can be divided in-place by instances of type `R` using
   * the `/=` operator.
   */
  template <typename L, typename R>
  concept DivideAssignAble = requires(L x, R y) { x /= y; };

  /**
   * \concept AddAssignAble
   * \brief Concept defining the requirements for types that support in-place addition.
   * \tparam L Type of the left operand.
   * \tparam R Type of the right operand.
   * \details The concept specifies that instances of type `L` can be added in-place by instances of type `R` using
   * the `+=` operator.
   */
  template <typename L, typename R>
  concept AddAssignAble = requires(L x, R y) { x += y; };

  /**
   * \concept SubstractAssignAble
   * \brief Concept defining the requirements for types that support in-place subtraction.
   * \tparam L Type of the left operand.
   * \tparam R Type of the right operand.
   * \details The concept specifies that instances of type `L` can be subtracted in-place by instances of type `R`
   * using the `-=` operator.
   */
  template <typename L, typename R>
  concept SubstractAssignAble = requires(L x, R y) { x -= y; };

  /**
   * \concept DivideAble
   * \brief Concept defining the requirements for types that support division.
   * \tparam L Type of the left operand.
   * \tparam R Type of the right operand.
   * \details The concept specifies that instances of types `L` and `R` can be divided using the `/` operator.
   */
  template <typename L, typename R>
  concept DivideAble = requires(L x, R y) { x / y; };

  /**
   * \concept NegateAble
   * \brief Concept defining the requirements for types that support negation.
   * \tparam L Type of the operand.
   * \details The concept specifies that instances of type `L` can be negated using the unary `-` operator.
   */
  template <typename L>
  concept NegateAble = requires(L x) { -x; };

  /**
   * \concept TransposeAble
   * \brief Concept defining the requirements for types that support transposition.
   * \tparam L Type of the operand.
   * \details The concept specifies that instances of type `L` can be transposed using the `transpose` function.
   */
  template <typename L>
  concept TransposeAble = requires(L x) { transpose(x); };

  /**
   * \concept IsFunctorWithArgs
   * \brief Concept defining the requirements for functors with arguments.
   * \tparam Op Type of the functor.
   * \tparam Args Types of the arguments.
   */
  template <typename Op, typename... Args>
  concept IsFunctorWithArgs = requires(Op op, Args... args) { op(args...); };

  /**
   * \concept EigenVector
   * \brief Concept defining the requirements for Eigen vectors.
   * \tparam V Type representing an Eigen vector.
   */
  template <typename V>
  concept EigenVector = static_cast<bool>(V::IsVectorAtCompileTime);

  /**
   * \concept EigenMatrix
   * \brief Concept defining the requirements for Eigen matrices. This also includes Eigen vectors
   * \tparam M Type representing an Eigen vector.
   */
  template <typename M>
  concept EigenMatrix = traits::isSpecializationTypeAndNonTypes<Eigen::Matrix, M>::value;

#define MAKE_EIGEN_FIXED_VECTOR_CONCEPT(Size) \
  template <typename V>                       \
  concept EigenVector##Size =                 \
      static_cast<bool>(V::IsVectorAtCompileTime) and static_cast<bool>(V::SizeAtCompileTime == Size);

  MAKE_EIGEN_FIXED_VECTOR_CONCEPT(1)
  MAKE_EIGEN_FIXED_VECTOR_CONCEPT(2)
  MAKE_EIGEN_FIXED_VECTOR_CONCEPT(3)
  MAKE_EIGEN_FIXED_VECTOR_CONCEPT(4)
  MAKE_EIGEN_FIXED_VECTOR_CONCEPT(5)
  MAKE_EIGEN_FIXED_VECTOR_CONCEPT(6)

#define MAKE_EIGEN_FIXED_MATRIX_CONCEPT(Size1, Size2)                                                           \
  template <typename M>                                                                                         \
  concept EigenMatrix##Size1##Size2 = static_cast<bool>(std::remove_cvref_t<M>::RowsAtCompileTime == Size1) and \
                                      static_cast<bool>(std::remove_cvref_t<M>::ColsAtCompileTime == Size2);

  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(1, 1)

  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(1, 2)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(2, 2)

  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(1, 3)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(2, 3)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(3, 2)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(3, 3)

  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(1, 4)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(2, 4)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(3, 4)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(4, 2)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(4, 3)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(4, 4)

  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(1, 5)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(2, 5)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(3, 5)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(4, 5)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(5, 2)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(5, 3)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(5, 4)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(5, 5)

  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(1, 6)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(2, 6)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(3, 6)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(4, 6)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(5, 6)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(6, 2)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(6, 3)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(6, 4)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(6, 5)
  MAKE_EIGEN_FIXED_MATRIX_CONCEPT(6, 6)

#define MAKE_EIGEN_FIXED_MATRIX_OR_VOIGT_CONCEPT(Size1, Size2) \
  template <typename M>                                        \
  concept EigenMatrixOrVoigtNotation##Size1 = EigenMatrix##Size1##Size1<M> or EigenVector##Size2<M>;

  MAKE_EIGEN_FIXED_MATRIX_OR_VOIGT_CONCEPT(1, 1)
  MAKE_EIGEN_FIXED_MATRIX_OR_VOIGT_CONCEPT(2, 3)
  MAKE_EIGEN_FIXED_MATRIX_OR_VOIGT_CONCEPT(3, 6)

  namespace Impl {
    template <template <typename...> class MaterialToCheck, typename Material>
    consteval bool isMaterial() {
      if constexpr (traits::isSpecialization<MaterialToCheck, Material>::value)
        return true;

      if constexpr (Material::isReduced) {
        if constexpr (traits::isSpecialization<MaterialToCheck, typename Material::Underlying>::value) {
          return true;
        } else {
          return false;
        }
      } else
        return false;
    }
  } // namespace Impl

  /**
   * \concept IsMaterial
   * \brief Concept defining the requirements for a material type.
   * \tparam MaterialToCheck Template representing the material type to check.
   * \tparam Material Type representing the material to be checked.
   * \details A type satisfies the IsMaterial concept if it meets one of the following conditions:
   *
   * 1. The material is a specialization of the specified template `MaterialToCheck`.
   * 2. The material is a specialization of `VanishingStress` or `VanishingStrain` (i.e. reduced) with an underlying
   * type that is a specialization of the specified template `MaterialToCheck`.
   *
   */
  template <template <typename...> class MaterialToCheck, typename Material>
  concept IsMaterial = Impl::isMaterial<MaterialToCheck, Material>();

  namespace Impl {
    template <typename T>
    concept ResultType = requires(T t) {
      typename T::type;                             // The nested type 'type'
      typename T::Vectorizer;                       // The nested type 'Vectorizer'
      typename T::Matricizer;                       // The nested type 'Matricizer'
      { toString(t) } -> std::same_as<std::string>; // The toString function
    };
  } // namespace Impl

  /**
   * \concept GeometricallyLinearMaterial
   * \brief Concepts defining the requirements for a material to be geometrically linear
   * This is the case when the corresponding strainTag is linear.
   * \tparam MAT the material implementation
   */
  template <typename MAT>
  concept GeometricallyLinearMaterial = MAT::strainTag == StrainTags::linear;

  /**
   * \concept ResultType
   * \brief A concept to check if a template type satisfies the ResultType requirements.
   * \tparam RT A template type with parameters (typename, int, int).
   *            The first parameter is the data type, and the next two parameters are the grid dimension and the world
   * dimension. It checks for various instantiations to ensure they meet the ResultType concept. Specifically, it checks
   * for the nested types 'type', 'Vectorizer', 'Matricizer', and the presence of a toString function returning a
   * std::string.
   */
  template <template <typename, int, int> typename RT>
  concept ResultType =
      Impl::ResultType<RT<double, 1, 1>> or Impl::ResultType<RT<double, 1, 2>> or Impl::ResultType<RT<double, 1, 3>> or
      Impl::ResultType<RT<double, 2, 3>> or Impl::ResultType<RT<double, 3, 3>>;

  /**

   * \brief Concept representing the requirements for a FlatAssembler.
   * \concept FlatAssembler
   * A type T satisfies FlatAssembler if it provides the necessary member functions
   * and data types for assembling sparse matrices in a flat structure.
   */
  template <typename T>
  concept FlatAssembler = requires(T t, const typename T::FEConfiguration& req,
                                   typename T::AffordanceCollectionType affordance, DBCOption dbcOption) {
    { t.requirement() } -> std::convertible_to<typename T::FEConfiguration&>;
    { t.affordanceCollection() } -> std::convertible_to<typename T::AffordanceCollectionType>;
    { t.dBCOption() } -> std::convertible_to<DBCOption>;

    { t.bind(req, affordance, dbcOption) } -> std::same_as<void>;
    { t.bind(req) } -> std::same_as<void>;
    { t.bind(affordance) } -> std::same_as<void>;
    { t.bind(dbcOption) } -> std::same_as<void>;

    { t.bound() } -> std::convertible_to<bool>;
    { t.boundToRequirement() } -> std::convertible_to<bool>;
    { t.boundToAffordanceCollection() } -> std::convertible_to<bool>;
    { t.boundToDBCOption() } -> std::convertible_to<bool>;
    { t.estimateOfConnectivity() } -> std::convertible_to<size_t>;

    { t.createFullVector(std::declval<Eigen::Ref<const Eigen::VectorXd>>()) } -> std::convertible_to<Eigen::VectorXd>;
    { t.constraintsBelow(std::declval<size_t>()) } -> std::convertible_to<size_t>;
    { t.isConstrained(std::declval<size_t>()) } -> std::convertible_to<bool>;
    { t.size() } -> std::convertible_to<size_t>;
    { t.reducedSize() } -> std::convertible_to<size_t>;
  };

  /**
   * \brief Concept representing the requirements for a ScalarFlatAssembler.
   * \concept ScalarFlatAssembler
   * A type T satisfies ScalarFlatAssembler if it is a FlatAssembler and if it provides the necessary scalar() member
   * functions.
   */
  template <typename T>
  concept ScalarFlatAssembler =
      Concepts::FlatAssembler<T> and requires(T t, const typename T::FEConfiguration& req,
                                              typename T::AffordanceCollectionType affordance, DBCOption dbcOption) {
        { t.scalar(req, affordance.scalarAffordance()) } -> std::convertible_to<const double&>;
        { t.scalar() } -> std::convertible_to<const double&>;
      };

  /**
   * \brief Concept representing the requirements for a VectorFlatAssembler.
   * \concept VectorFlatAssembler
   * A type T satisfies VectorFlatAssembler if it is a ScalarFlatAssembler and if it provides the necessary vector()
   * member functions.
   */
  template <typename T>
  concept VectorFlatAssembler = Concepts::ScalarFlatAssembler<T> and
                                requires(T t, const typename T::FEConfiguration& req,
                                         typename T::AffordanceCollectionType affordance, DBCOption dbcOption) {
                                  {
                                    t.vector(req, affordance.vectorAffordance(), dbcOption)
                                  } -> std::convertible_to<const Eigen::VectorXd&>;
                                  { t.vector(dbcOption) } -> std::convertible_to<const Eigen::VectorXd&>;
                                  { t.vector() } -> std::convertible_to<const Eigen::VectorXd&>;
                                };

  /**
   * \brief Concept representing the requirements for a MatrixFlatAssembler.
   * \concept MatrixFlatAssembler
   * A type T satisfies MatrixFlatAssembler if it is a VectorFlatAssembler and if it provides the necessary matrix()
   * member functions.
   */
  template <typename T>
  concept MatrixFlatAssembler = Concepts::VectorFlatAssembler<T> and
                                requires(T t, const typename T::FEConfiguration& req,
                                         typename T::AffordanceCollectionType affordance, DBCOption dbcOption) {
                                  { t.matrix(req, affordance.matrixAffordance(), dbcOption) };
                                  { t.matrix(dbcOption) };
                                  { t.matrix() };
                                };

  // adapted from /dune/dune-vtk/dune/vtk/utility/concepts.hh
  template <class DC>
  concept DataCollector = requires(DC dc) {
    typename DC::GridView;
    { dc.update() } -> std::same_as<void>;
    { dc.numPoints() } -> std::convertible_to<std::uint64_t>;
    { dc.numCells() } -> std::convertible_to<std::uint64_t>;
    { dc.gridView() } -> std::same_as<const typename DC::GridView&>;
  };

  template <class GV>
  concept GridView = requires(GV g) {
    typename GV::Grid;
    GV::dimension;
    GV::dimensionworld;

    { g.grid() };
  };

  namespace Impl {
    template <typename T>
    struct is_dual : std::false_type
    {
    };

    // Specialization for Dual<T, U>: this will be true for Dual types
    template <typename T, typename U>
    struct is_dual<autodiff::detail::Dual<T, U>> : std::true_type
    {
    };
  } // namespace Impl

  /**
  \concept AutodiffScalar
  \brief Concept to check if the underlying scalar type is a dual type.
  \tparam T The scalar type to be checked if it is a dual type.
  */
  template <typename T>
  concept AutodiffScalar = Impl::is_dual<T>::value;

} // namespace Concepts
} // namespace Ikarus
