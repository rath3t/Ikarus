// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <type_traits>

#include <Eigen/Core>

#include <ikarus/utils/tensorutils.hh>

/**
 * \file feresulttypes.hh
 * \brief Definitions of ResultTypes used for finite element results
 * \ingroup finiteelements
 * \details Additional ResultTypes can be added with one of the following macros
 */
namespace Ikarus {
namespace Impl {

  template <bool strainlike = false>
  struct VectorizeWithVoigt
  {
    template <typename Derived>
    static auto transform(const Eigen::DenseBase<Derived>& mat) {
      return toVoigt(mat.derived(), strainlike);
    }
  };

  struct VectorizeGeneric
  {
    template <typename Derived>
    static auto transform(const Eigen::DenseBase<Derived>& mat) {
      return mat.derived().reshaped().eval();
    }
  };

  template <bool strainlike = false>
  struct MatricizeWithVoigt
  {
    template <typename Derived, int RowsAtCompileTime, int ColsAtCompileTime>
    static auto transform(const Eigen::DenseBase<Derived>& vec, int rows = RowsAtCompileTime,
                          int cols = ColsAtCompileTime) {
      assert(rows == RowsAtCompileTime && cols == ColsAtCompileTime &&
             "Only the fixed size values work for voigt matrices and vectors");
      static_assert(RowsAtCompileTime != Eigen::Dynamic and ColsAtCompileTime != Eigen::Dynamic,
                    "Voigt notation only available for fixed size vectors and matrices");
      return fromVoigt(vec.derived(), strainlike);
    }
  };

  struct MatricizeGeneric
  {
    template <typename Derived, int RowsAtCompileTime, int ColsAtCompileTime>
    static auto transform(const Eigen::DenseBase<Derived>& vec, int rows = RowsAtCompileTime,
                          int cols = ColsAtCompileTime) {
      return vec.derived().reshaped(Eigen::fix<RowsAtCompileTime>(rows), Eigen::fix<ColsAtCompileTime>(cols)).eval();
    }
  };

} // namespace Impl

namespace ResultTypes {
#define REGISTER_RESULTTYPE_IMPL(resultTypeName, rowsExpr, colsExpr, MaxRowsExpr, MaxColsExpr, VectorizeStruct, \
                                 MatricizeStruct)                                                               \
  template <typename ScalarType, int gridDim, int worldDim>                                                     \
  struct resultTypeName                                                                                         \
  {                                                                                                             \
    friend std::string toString(resultTypeName) { return #resultTypeName; }                                     \
                                                                                                                \
    using type       = Eigen::Matrix<ScalarType, rowsExpr, colsExpr, 0, MaxRowsExpr, MaxColsExpr>;              \
    using Vectorizer = VectorizeStruct;                                                                         \
    using Matricizer = MatricizeStruct;                                                                         \
                                                                                                                \
    template <typename ScalarType_, int gridDim_, int worldDim_>                                                \
    using Rebind = resultTypeName<ScalarType_, gridDim_, worldDim_>;                                            \
  }

/**
 * \brief Used to register a symmetric ResultType with compile-time fixed rows and columns (uses Voigt notation)
 * \param resultTypeName name of the ResultType
 * \param rowsExpr expression for rows, e.g. `gridDim` or `worldDim`
 * \param colsExpr expression for rows, e.g. `gridDim` or `worldDim`
 * \param strainlike boolean indicating whether the result should be treated as a strain-like quantity
 */
#define REGISTER_SIMPLE_SYMMETRIC_RESULTTYPE(resultTypeName, rowsExpr, colsExpr, strainlike) \
  REGISTER_RESULTTYPE_IMPL(resultTypeName, rowsExpr, colsExpr, rowsExpr, colsExpr,           \
                           Ikarus::Impl::VectorizeWithVoigt<strainlike>, Ikarus::Impl::MatricizeWithVoigt<strainlike>)

  /**
   * \brief Used to register a general ResultType with potentially dynamic size without reserved memeroy
   * rows and columns
   * \param resultTypeName name of the ResultType
   * \param rowsExpr expression for rows, e.g. `Eigen::Dynamic` or `worldDim`
   * \param colsExpr expression for columns, e.g. `Eigen::Dynamic` or `worldDim`
   */
#define REGISTER_RESULTTYPE(resultTypeName, rowsExpr, colsExpr)                                \
  REGISTER_RESULTTYPE_IMPL(resultTypeName, rowsExpr, colsExpr, Ikarus::Impl::VectorizeGeneric, \
                           Ikarus::Impl::MatricizeGeneric)
  /**
   * \brief Used to register a general ResultType with potentially dynamic size and defined maximum amount of
   * rows and columns
   * \param resultTypeName name of the ResultType
   * \param rowsExpr expression for rows, e.g. `Eigen::Dynamic` or `worldDim`
   * \param colsExpr expression for columns, e.g. `Eigen::Dynamic` or `worldDim`
   * \param MaxRowsExpr expression for maximum number of columns
   * \param MaxColsExprexpression for maximum number of rows
   */
#define REGISTER_RESERVED_RESULTTYPE(resultTypeName, rowsExpr, colsExpr, MaxRowsExpr, MaxColsExpr) \
  REGISTER_RESULTTYPE_IMPL(resultTypeName, rowsExpr, colsExpr, MaxRowsExpr, MaxColsExpr,           \
                           Ikarus::Impl::VectorizeGeneric, Ikarus::Impl::MatricizeGeneric)
/**
 * \brief Used to register a general ResultType with potentially dynamic size, with reserved memory according to
 * passed in rowsExpr and colsExpr
 * \param resultTypeName name of the ResultType
 * \param rowsExpr expression for rows, e.g. `Eigen::Dynamic` or `worldDim`
 * \param colsExpr expression for columns, e.g. `Eigen::Dynamic` or `worldDim`
 */
#define REGISTER_SIMPLE_RESULTTYPE(resultTypeName, rowsExpr, colsExpr) \
  REGISTER_RESERVED_RESULTTYPE(resultTypeName, rowsExpr, colsExpr, rowsExpr, colsExpr)

  ///=============================
  /// Registered ResultTypes

  REGISTER_SIMPLE_SYMMETRIC_RESULTTYPE(linearStress, worldDim, worldDim, false);
  REGISTER_SIMPLE_SYMMETRIC_RESULTTYPE(PK2Stress, worldDim, worldDim, false);
  REGISTER_SIMPLE_SYMMETRIC_RESULTTYPE(cauchyStress, worldDim, worldDim, false);

  REGISTER_SIMPLE_RESULTTYPE(director, worldDim, 1);
  REGISTER_SIMPLE_RESULTTYPE(magnetization, worldDim, 1);
  REGISTER_SIMPLE_RESULTTYPE(gradientNormOfMagnetization, 1, 1);
  REGISTER_SIMPLE_RESULTTYPE(vectorPotential, worldDim, 1);
  REGISTER_SIMPLE_RESULTTYPE(divergenceOfVectorPotential, 1, 1);

  REGISTER_SIMPLE_RESULTTYPE(BField, worldDim, 1);
  REGISTER_SIMPLE_RESULTTYPE(HField, worldDim, 1);

  REGISTER_SIMPLE_RESULTTYPE(cauchyAxialForce, 1, 1);
  REGISTER_SIMPLE_RESULTTYPE(PK2AxialForce, 1, 1);
  REGISTER_SIMPLE_RESULTTYPE(linearAxialForce, 1, 1);

  REGISTER_SIMPLE_RESULTTYPE(customType, Eigen::Dynamic, Eigen::Dynamic);
} // namespace ResultTypes

enum class ResultShape
{
  Vector,
  Matrix
};

/**
 * \brief Container that is used for FE Results. It gives access to the stored value, but can also be used to access
 * the result in Matrix and Vector form
 * \tparam RT A specified ResultType
 * \tparam storedResultShape boolean indicating whether the stored result is in its vector form, defaults to true
 */
template <typename RT, ResultShape storedResultShape = ResultShape::Vector>
struct ResultWrapper : RT
{
private:
  using ResultTypeValueType                       = typename RT::type;
  static constexpr Eigen::Index rowsAtCompileTime = ResultTypeValueType::RowsAtCompileTime;
  static constexpr Eigen::Index colsAtCompileTime = ResultTypeValueType::ColsAtCompileTime;
  static constexpr bool storedValueIsVector       = storedResultShape == ResultShape::Vector;

public:
  using VecType = std::invoke_result_t<decltype(&RT::Vectorizer::template transform<ResultTypeValueType>),
                                       const ResultTypeValueType&>;
  using MatType =
      std::invoke_result_t<decltype(&RT::Matricizer::template transform<VecType, rowsAtCompileTime, colsAtCompileTime>),
                           const VecType&, int, int>;
  using StoredType = std::conditional_t<storedValueIsVector, VecType, MatType>;
  using ResultType = RT;

  /**
   * \brief Returns the stored value as Vector
   * \return vector representation of the stored value
   */
  auto asVec() const {
    if constexpr (storedValueIsVector)
      return value_;
    else
      return RT::Vectorizer::transform(value_);
  }

  /**
   * \brief Returns the stored value as Matrix (if possible)
   * \param rows specify the rows of the matrix for dynamic sized Result Types
   * \param cols specify the columns of the matrix for dynamic sized Result Types
   * \return matrix representation of the stored value
   */
  auto asMat(Eigen::Index rows = rowsAtCompileTime, Eigen::Index cols = colsAtCompileTime) const {
    if constexpr (storedValueIsVector) {
      if constexpr (rowsAtCompileTime == Eigen::Dynamic)
        assert(rows != rowsAtCompileTime &&
               "For dynamic size result types you have to pass rows by hand, since it is not clear how the result "
               "should be reshaped");
      if constexpr (colsAtCompileTime == Eigen::Dynamic)
        assert(cols != colsAtCompileTime &&
               "For dynamic size result types you have to pass cols by hand, since it is not clear how the result "
               "should be reshaped");
      return RT::Matricizer::template transform<VecType, rowsAtCompileTime, colsAtCompileTime>(value_, rows, cols);
    } else
      return value_;
  }
  explicit ResultWrapper() = default;
  explicit ResultWrapper(StoredType&& value) { this->value_ = std::move(value); }
  explicit ResultWrapper(const StoredType& value) { this->value_ = value; }
  ResultWrapper& operator=(const StoredType& value) {
    this->value_ = value;
    return *this;
  }
  ResultWrapper& operator=(StoredType&& value) {
    this->value_ = std::move(value);
    return *this;
  }

private:
  StoredType value_{};
};

namespace Impl {
  template <template <typename, int, int> class RT>
  using DummyRT = RT<double, 1, 1>;
}

/**
 * \brief Creates a dummy resultType which can be stored in a variable
 * \tparam RT The ResultType template
 * \return The ResultType
 */
template <template <typename, int, int> class RT>
auto makeRT() {
  return Impl::DummyRT<RT>{};
}

/**
 * \brief Retrieves a string representation of the ResultType template
 * \tparam RT the ResultType template
 * \return the name of the ResultType
 */
template <template <typename, int, int> class RT>
auto toString() {
  return toString(Impl::DummyRT<RT>{});
}

/**
 * \brief Meta variable to test whether two ResultType templates are the same
 * \tparam RT1 first ResultType template
 * \tparam RT2 second ResultType template
 */
template <template <typename, int, int> class RT1, template <typename, int, int> class RT2>
constexpr bool isSameResultType = std::is_same_v<Impl::DummyRT<RT1>, Impl::DummyRT<RT2>>;

namespace Impl {
  template <typename T, typename Tuple>
  struct hasType;

  template <typename T, typename... Us>
  struct hasType<T, std::tuple<Us...>> : std::disjunction<std::is_same<T, Us>...>
  {
  };
} // namespace Impl

/**
 * \brief Base class for element definitions that provides common functionality for ResultTypes
 *
 * \tparam ResultTypes supported ResultTypes
 */
template <template <typename, int, int> typename... ResultTypes>
struct ResultTypeBase
{
  /**
   * \brief Returns whether a ResultType is provided by the element
   * \tparam RT requested ResultType
   */
  template <template <typename, int, int> typename RT>
  static consteval bool supportsResultType() {
    return Impl::hasType<decltype(makeRT<RT>()), SupportedResultTypes>::value;
  }

  using SupportedResultTypes = std::tuple<decltype(makeRT<ResultTypes>())...>;
};

} // namespace Ikarus
