// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file linearalgebrahelper.hh
 * \brief Helper for the autodiff library
 */

#pragma once
#include "concepts.hh"
#include "traits.hh"

#include <iosfwd>
#include <random>
#include <ranges>

#include <dune/common/tuplevector.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/istl/bvector.hh>

#include <Eigen/Core>

#include <autodiff/forward/dual/dual.hpp>

namespace Ikarus {

/**
 * \brief Orthonormalizes all Matrix columns using Gram-Schmidt Orthogonalization.
 * \tparam Derived Type of the input matrix.
 * \ingroup utils
 * \param A The input matrix.
 * \return Eigen Matrix with orthonormalized columns.
 */
template <typename Derived>
auto orthonormalizeMatrixColumns(const Eigen::MatrixBase<Derived>& A) {
  // Gram Schmidt Ortho
  auto Q = A.eval();

  Q.col(0).normalize();

  for (int colIndex = 1; colIndex < Q.cols(); colIndex++) {
    Q.col(colIndex) -= Q.leftCols(colIndex) * (Q.leftCols(colIndex).transpose() * A.col(colIndex));
    Q.col(colIndex).normalize();
  }

  return Q;
}

/**
 * \brief View Dune::BlockVector as an Eigen::Vector.
 * \ingroup utils
 * \tparam ValueType Type of the elements in the BlockVector.
 * \param blockedVector Input Dune::BlockVector.
 * \return Eigen::Map of the BlockVector as a flat Eigen::Vector.
 */
template <typename ValueType>
auto viewAsFlatEigenVector(Dune::BlockVector<ValueType>& blockedVector) {
  Eigen::Map<Eigen::VectorX<typename ValueType::field_type>> vec(&blockedVector.begin()->begin().operator*(),
                                                                 blockedVector.size() * blockedVector[0].size());

  return vec;
}

/**
 * \brief View const Dune::BlockVector as an Eigen::Vector.
 *  \ingroup utils
 * \tparam ValueType Type of the elements in the BlockVector.
 * \param blockedVector Input Dune::BlockVector.
 * \return Eigen::Map of the BlockVector as a flat Eigen::Vector (const version).
 */
template <typename ValueType>
auto viewAsFlatEigenVector(const Dune::BlockVector<ValueType>& blockedVector) {
  Eigen::Map<const Eigen::VectorX<typename ValueType::field_type>> vec(&blockedVector.begin()->begin().operator*(),
                                                                       blockedVector.size() * blockedVector[0].size());

  return vec;
}

/**
 * \brief View Dune::BlockVector as an Eigen::Matrix with dynamic rows and fixed columns depending on the size of the
 * ValueType.
 * \ingroup utils
 * \tparam ValueType Type of the elements in the BlockVector.
 * \param blockedVector Input Dune::BlockVector.
 * \return Eigen::Map of the BlockVector as a dynamic Eigen::Matrix.
 */
template <typename ValueType>
auto viewAsEigenMatrixAsDynFixed(Dune::BlockVector<ValueType>& blockedVector) {
  Eigen::Map<Eigen::Matrix<typename ValueType::field_type, Eigen::Dynamic, ValueType::valueSize, Eigen::RowMajor>> vec(
      &blockedVector.begin()->begin().operator*(), blockedVector.size(), blockedVector[0].size());

  return vec;
}

/**
 * \brief Const view Dune::BlockVector as an Eigen::Matrix with dynamic rows and fixed columns depending on the size
 * of the ValueType.
 * \ingroup utils
 * \tparam ValueType Type of the elements in the BlockVector.
 * \param blockedVector Input Dune::BlockVector.
 * \return Eigen::Map of the BlockVector as a dynamic Eigen::Matrix (const version).
 */
template <typename ValueType>
auto viewAsEigenMatrixAsDynFixed(const Dune::BlockVector<ValueType>& blockedVector) {
  Eigen::Map<const Eigen::Matrix<typename ValueType::field_type, Eigen::Dynamic, ValueType::valueSize, Eigen::RowMajor>>
      vec(&blockedVector.begin()->begin().operator*(), blockedVector.size(), blockedVector[0].size());

  return vec;
}

/**
 * \brief View Dune::BlockVector as an Eigen::Matrix with fixed rows depending on the size of the ValueType and
 * dynamic columns.
 * \ingroup utils
 * \tparam ValueType Type of the elements in the BlockVector.
 * \param blockedVector Input Dune::BlockVector.
 * \return Eigen::Map of the BlockVector as a fixed-size Eigen::Matrix with dynamic columns.
 */
template <typename ValueType>
auto viewAsEigenMatrixFixedDyn(Dune::BlockVector<ValueType>& blockedVector) {
  Eigen::Map<Eigen::Matrix<typename ValueType::field_type, ValueType::valueSize, Eigen::Dynamic>> vec(
      &blockedVector.begin()->begin().operator*(), blockedVector[0].size(), blockedVector.size());

  return vec;
}

/**
 * \brief Const view Dune::BlockVector as an Eigen::Matrix with fixed rows depending on the size of the ValueType and
 * dynamic columns. \ingroup utils
 * \tparam ValueType Type of the elements in the BlockVector.
 * \param blockedVector Input Dune::BlockVector.
 * \return Eigen::Map of the BlockVector as a fixed-size Eigen::Matrix with dynamic columns (const version).
 */
template <typename ValueType>
auto viewAsEigenMatrixFixedDyn(const Dune::BlockVector<ValueType>& blockedVector) {
  Eigen::Map<const Eigen::Matrix<typename ValueType::field_type, ValueType::valueSize, Eigen::Dynamic>> vec(
      &blockedVector.begin()->begin().operator*(), blockedVector[0].size(), blockedVector.size());

  return vec;
}

/**
 * \brief Returns the total correction size of a block vector with a Manifold as the underlying type.
 *  \ingroup utils
 * \tparam Type Manifold type.
 * \param a Input Dune::BlockVector.
 * \return Total correction size.
 */
template <typename Type>
size_t correctionSize(const Dune::BlockVector<Type>& a)
requires requires { Type::correctionSize; }
{
  return a.size() * Type::correctionSize;
}

/**
 * \brief Returns the total value size of a block vector with a Manifold as the underlying type.
 *  \ingroup utils
 * \tparam Type Manifold type.
 * \param a Input Dune::BlockVector.
 * \return Total value size.
 */
template <typename Type>
size_t valueSize(const Dune::BlockVector<Type>& a)
requires requires { Type::valueSize; }
{
  return a.size() * Type::valueSize;
}

/**
 * \brief Enables the += operator for Dune::BlockVector += Eigen::Vector.
 *  \ingroup utils
 * \tparam Type Manifold type.
 * \tparam Derived Type of the input Eigen matrix.
 * \param a Input Dune::BlockVector.
 * \param b Input Eigen::Matrix.
 * \return Reference to the modified Dune::BlockVector.
 */
template <typename Type, typename Derived>
Dune::BlockVector<Type>& operator+=(Dune::BlockVector<Type>& a, const Eigen::MatrixBase<Derived>& b)
requires(Ikarus::Concepts::AddAssignAble<Type, decltype(b.template segment<Type::correctionSize>(0))> and
         requires() { Type::correctionSize; })
{
  assert(correctionSize(a) == static_cast<size_t>(b.size()) && " The passed vector has wrong size");
  for (auto i = 0U; i < a.size(); ++i)
    a[i] += b.template segment<Type::correctionSize>(i * Type::correctionSize);
  return a;
}

/**
 * \brief Enables the -= operator for Dune::BlockVector += Eigen::Vector.
 *  \ingroup utils
 * \tparam Type Manifold type.
 * \tparam Derived Type of the input Eigen matrix.
 * \param a Input Dune::BlockVector.
 * \param b Input Eigen::Matrix.
 * \return Reference to the modified Dune::BlockVector.
 */
template <typename Type, typename Derived>
Dune::BlockVector<Type>& operator-=(Dune::BlockVector<Type>& a, const Eigen::MatrixBase<Derived>& b)
requires(Ikarus::Concepts::AddAssignAble<Type, decltype(b.template segment<Type::correctionSize>(0))> and
         requires() { Type::correctionSize; })
{
  return a += (-b);
}

/**
 * \brief Enables the += operator for Dune::TupleVector += Eigen::Vector.
 * \ingroup utils
 * \tparam Types Types of the elements in the TupleVector.
 * \tparam Derived Type of the input Eigen matrix.
 * \param a Input Dune::TupleVector.
 * \param b Input Eigen::Matrix.
 * \return Reference to the modified Dune::TupleVector.
 */
template <typename... Types, typename Derived>
Dune::TupleVector<Types...>& operator+=(Dune::TupleVector<Types...>& a, const Eigen::MatrixBase<Derived>& b) {
  using namespace Dune::Indices;
  size_t posStart = 0;
  Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<sizeof...(Types)>()), [&](const auto i) {
    const size_t size = correctionSize(a[i]);
    a[i] += b(Eigen::seqN(posStart, size));
    posStart += size;
  });

  return a;
}

/**
 * \brief Enables the addition in the embedding space of a vector in the space M^n, where M is a manifold with the
 * points of type ManifoldPoint \ingroup utils
 * \tparam ManifoldPoint Manifold type.
 * \tparam Derived ManifoldPoint of the input Eigen matrix.
 * \param a Input Dune::BlockVector.
 * \param b Input Eigen::Matrix.
 * \return Reference to the modified Dune::BlockVector.
 */
template <typename ManifoldPoint, typename Derived>
Dune::BlockVector<ManifoldPoint>& addInEmbedding(Dune::BlockVector<ManifoldPoint>& a,
                                                 const Eigen::MatrixBase<Derived>& b)
requires(Ikarus::Concepts::AddAssignAble<ManifoldPoint, decltype(b.template segment<ManifoldPoint::valueSize>(0))> and
         requires() { ManifoldPoint::valueSize; })
{
  assert(valueSize(a) == static_cast<size_t>(b.size()) && " The passed vector has wrong size");
  for (auto i = 0U; i < a.size(); ++i)
    a[i].addInEmbedding(b.template segment<ManifoldPoint::valueSize>(i * ManifoldPoint::valueSize));
  return a;
}

/**
 * \brief Adding free norm function to Eigen types.
 *  \ingroup utils
 * \tparam Derived Type of the input Eigen matrix.
 * \param v Input Eigen matrix.
 * \return Norm of the matrix.
 */
template <typename Derived>
requires(!std::floating_point<Derived>)
auto norm(const Eigen::MatrixBase<Derived>& v) {
  return v.norm();
}

/**
 * \brief Helper Free Function to have the same interface as for Eigen Vector Types.
 *  \ingroup utils
 * \param v Input scalar.
 * \return Absolute value of the scalar.
 */
auto norm(const std::floating_point auto& v) { return std::abs(v); }

/**
 * \brief Eigen::DiagonalMatrix Product Missing in Eigen.
 *  \ingroup utils
 * \tparam Scalar Scalar type.
 * \tparam size Size of the diagonal matrix.
 * \param a Input DiagonalMatrix.
 * \param b Input DiagonalMatrix.
 * \return Product of the two DiagonalMatrices.
 */
template <typename Scalar, int size>
auto operator*(const Eigen::DiagonalMatrix<Scalar, size>& a, const Eigen::DiagonalMatrix<Scalar, size>& b) {
  return (a.diagonal().cwiseProduct(b.diagonal())).asDiagonal();
}

/**
 * \brief In-place addition for Eigen::DiagonalMatrix.
 *  \ingroup utils
 * \tparam Scalar Scalar type.
 * \tparam size Size of the diagonal matrix.
 * \param a Input DiagonalMatrix.
 * \param b Input DiagonalMatrix.
 * \return Reference to the modified DiagonalMatrix.
 */
template <typename Scalar, int size>
auto operator+=(Eigen::DiagonalMatrix<Scalar, size>& a, const Eigen::DiagonalMatrix<Scalar, size>& b) {
  a.diagonal() += b.diagonal();
  return a;
}

/**
 * \brief Eigen::Matrix + Eigen::DiagonalMatrix addition missing in Eigen.
 *  \ingroup utils
 * \tparam Derived Type of the input Eigen matrix.
 * \tparam Scalar Scalar type.
 * \tparam size Size of the diagonal matrix.
 * \param a Input Eigen matrix.
 * \param b Input DiagonalMatrix.
 * \return Sum of the Eigen matrix and DiagonalMatrix.
 */
template <typename Derived, typename Scalar, int size>
auto operator+(const Eigen::MatrixBase<Derived>& a, const Eigen::DiagonalMatrix<Scalar, size>& b) {
  auto c = a.derived().eval();
  c.diagonal() += b.diagonal();
  return c;
}

/**
 * \brief Eigen::DiagonalMatrix + Eigen::Matrix addition missing in Eigen.
 *  \ingroup utils
 * \tparam Derived Type of the input Eigen matrix.
 * \tparam Scalar Scalar type.
 * \tparam size Size of the diagonal matrix.
 * \param a Input DiagonalMatrix.
 * \param b Input Eigen matrix.
 * \return Sum of the DiagonalMatrix and Eigen matrix.
 */
template <typename Derived, typename Scalar, int size>
auto operator+(const Eigen::DiagonalMatrix<Scalar, size>& a, const Eigen::MatrixBase<Derived>& b) {
  return b + a;
}

/**
 * \brief Unary minus for Eigen::DiagonalMatrix.
 *  \ingroup utils
 * \tparam Scalar Scalar type.
 * \tparam size Size of the diagonal matrix.
 * \param a Input DiagonalMatrix.
 * \return Negation of the DiagonalMatrix.
 */
template <typename Scalar, int size>
auto operator-(const Eigen::DiagonalMatrix<Scalar, size>& a) {
  return (-a.diagonal()).asDiagonal();
}

/**
 * \brief Addition of Eigen::Matrix and Eigen::DiagonalWrapper.
 *  \ingroup utils
 * \tparam Derived Type of the input Eigen matrix.
 * \tparam Derived2 Type of the input Eigen DiagonalWrapper.
 * \param a Input Eigen matrix.
 * \param b Input Eigen DiagonalWrapper.
 * \return Sum of Eigen matrix and DiagonalWrapper.
 */
template <typename Derived, typename Derived2>
auto operator+(const Eigen::MatrixBase<Derived>& a, const Eigen::DiagonalWrapper<Derived2>& b) {
  auto c = a.derived().eval();
  c.diagonal() += b.diagonal();
  return c;
}

/**
 * \brief Addition of Eigen::DiagonalWrapper and Eigen::Matrix.
 *  \ingroup utils
 * \tparam Derived Type of the input Eigen DiagonalWrapper.
 * \tparam Derived2 Type of the input Eigen matrix.
 * \param a Input Eigen DiagonalWrapper.
 * \param b Input Eigen matrix.
 * \return Sum of DiagonalWrapper and Eigen matrix.
 */
template <typename Derived, typename Derived2>
auto operator+(const Eigen::DiagonalWrapper<Derived>& a, const Eigen::MatrixBase<Derived2>& b) {
  return b + a;
}

/**
 * \brief Output stream operator for Eigen::DiagonalMatrix.
 *  \ingroup utils
 * \tparam Scalar Scalar type.
 * \tparam size Size of the diagonal matrix.
 * \param os Output stream.
 * \param a Input DiagonalMatrix.
 * \return Reference to the output stream.
 */
template <typename Scalar, int size>
std::ostream& operator<<(std::ostream& os, const Eigen::DiagonalMatrix<Scalar, size>& a) {
  os << Eigen::Matrix<Scalar, size, size>(a);
  return os;
}

/**
 * \brief Returns the symmetric part of a matrix.
 *  \ingroup utils
 * \tparam Derived Type of the input Eigen matrix.
 * \param A Input Eigen matrix.
 * \return Symmetric part of the matrix.
 */
template <typename Derived>
Derived sym(const Eigen::MatrixBase<Derived>& A) {
  return 0.5 * (A + A.transpose());
}

/**
 * \brief Returns the skew part of a matrix.
 *  \ingroup utils
 * \tparam Derived Type of the input Eigen matrix.
 * \param A Input Eigen matrix.
 * \return Skew part of the matrix.
 */
template <typename Derived>
Derived skew(const Eigen::MatrixBase<Derived>& A) {
  return 0.5 * (A - A.transpose());
}

/**
 * \brief Method to print the matrix in a format that can directly be copied to Maple.
 *  \ingroup utils
 * \tparam Derived The derived type of the matrix.
 * \param A The input matrix.
 */
template <typename Derived>
void printForMaple(const Eigen::EigenBase<Derived>& A) {
  Eigen::IOFormat mapleFmt(Eigen::FullPrecision, 0, ", ", "|\n", "<", ">", "<", ">");
  if constexpr (std::convertible_to<Derived, const Eigen::MatrixBase<Derived>&>) {
    std::cout << "\n" << A.derived().format(mapleFmt) << std::endl;
  } else { // branch for Eigen::DiagonalMatrix
    using Scalar = typename Derived::Scalar;
    using namespace Eigen;
    constexpr int diag_size = EIGEN_SIZE_MIN_PREFER_DYNAMIC(Derived::RowsAtCompileTime, Derived::ColsAtCompileTime);
    std::cout << "\n"
              << Eigen::Matrix<Scalar, diag_size, diag_size>(A.derived().diagonal().asDiagonal()).format(mapleFmt)
              << std::endl;
  }
}

/**
 * \brief Creates a random vector of the specified type within a given range.
 *  \ingroup utils
 * \tparam FieldVectorT The type of the vector.
 * \param lower The lower bound of the random values (default is -1).
 * \param upper The upper bound of the random values (default is 1).
 * \return A random vector within the specified range.
 */
template <typename FieldVectorT>
auto createRandomVector(typename FieldVectorT::value_type lower = -1, typename FieldVectorT::value_type upper = 1) {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<typename FieldVectorT::value_type> dist(lower, upper);
  auto rand = [&dist, &mt]() { return dist(mt); };
  FieldVectorT vec;
  std::generate(vec.begin(), vec.end(), rand);
  return vec;
}

/**
 * \brief Create skew 3x3 matrix from 3d vector.
 *  \ingroup utils
 * \tparam ScalarType The type of the coordinates in the vector.
 * \param a The vector.
 * \return The skew matrix.
 */
template <typename ScalarType>
Eigen::Matrix<ScalarType, 3, 3> skew(const Eigen::Vector<ScalarType, 3>& a) {
  Eigen::Matrix<ScalarType, 3, 3> A;
  A << 0, -a(2), a(1), a(2), 0, -a(0), -a(1), a(0), 0;
  return A;
}

namespace Impl {
  constexpr std::tuple<std::array<std::array<int, 2>, 1>, std::array<std::array<int, 2>, 3>,
                       std::array<std::array<int, 2>, 6>>
      voigtIndices = {{{{0, 0}}}, {{{0, 0}, {1, 1}, {0, 1}}}, {{{0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}}}};
}

/**
 * \brief Container for Voigt notation indices based on dimension.
 *  \ingroup utils
 * 1D: 0,0
 * 2D: 0,0; 1,1; 0,1
 * 3D: 0,0; 1,1; 2,2; 1,2; 0,2; 0,1
 * \tparam dim The dimension for which Voigt indices are needed.
 */
template <int dim>
constexpr auto voigtNotationContainer = std::get<dim - 1>(Impl::voigtIndices);

/**
 * \brief Performs static condensation on a square matrix.
 *  \ingroup utils
 * \tparam Derived Type of the matrix.
 * \tparam sizeOfCondensedIndices Size of the condensed indices.
 * \param E Input matrix.
 * \param indices Array of indices to be condensed.
 * \return Resulting matrix after static condensation.
 *  \details
 * This function performs static condensation on a square matrix. It removes the specified indices from the matrix,
 * computes the remaining submatrices (K11, K12, K22), and returns the result of the static condensation.
 */
template <typename Derived, size_t sizeOfCondensedIndices>
auto staticCondensation(const Eigen::MatrixBase<Derived>& E,
                        const std::array<size_t, sizeOfCondensedIndices>& indices) {
  constexpr size_t colsFull = std::remove_cvref_t<Derived>::ColsAtCompileTime;
  constexpr size_t rowsFull = std::remove_cvref_t<Derived>::RowsAtCompileTime;
  static_assert(colsFull == rowsFull, "staticCondensation only allowed for square matrices");
  std::array<size_t, rowsFull - sizeOfCondensedIndices> remainingIndices{};
  std::ranges::set_difference(std::ranges::iota_view(size_t(0), size_t(colsFull)), indices, remainingIndices.begin());

  auto K11 = E(remainingIndices, remainingIndices);
  auto K12 = E(indices, remainingIndices);
  auto K22 = E(indices, indices);

  return (K11 - K12.transpose() * K22.inverse() * K12).eval();
}

template <typename Derived, size_t sizeOfCondensedIndices>
auto reduceMatrix(const Eigen::MatrixBase<Derived>& E, const std::array<size_t, sizeOfCondensedIndices>& indices) {
  constexpr size_t colsFull = std::remove_cvref_t<Derived>::ColsAtCompileTime;
  constexpr size_t rowsFull = std::remove_cvref_t<Derived>::RowsAtCompileTime;
  static_assert(colsFull == rowsFull, "reduceMatrix only allowed for square matrices");
  std::array<size_t, rowsFull - sizeOfCondensedIndices> remainingIndices{};
  std::ranges::set_difference(std::ranges::iota_view(size_t(0), size_t(colsFull)), indices, remainingIndices.begin());

  auto K11 = E(remainingIndices, remainingIndices);

  return K11.eval();
}

/**
 * \brief Removes specified columns from a matrix.
 *  \ingroup utils
 * \tparam Derived Type of the matrix.
 * \tparam sizeOfRemovedCols Size of the columns to be removed.
 * \param E Input matrix.
 * \param indices Array of column indices to be removed.
 * \return Resulting matrix after removing specified columns.
 *  \details
 * This function removes specified columns from a matrix. It computes the remaining columns after removing the
 * specified indices and returns the resulting matrix.
 */
template <typename Derived, size_t sizeOfRemovedCols>
auto removeCol(const Eigen::MatrixBase<Derived>& E, const std::array<size_t, sizeOfRemovedCols>& indices) {
  constexpr size_t colsFull = std::remove_cvref_t<Derived>::ColsAtCompileTime;
  constexpr size_t rowsFull = std::remove_cvref_t<Derived>::RowsAtCompileTime;
  static_assert(colsFull == 1);

  std::array<size_t, rowsFull - sizeOfRemovedCols> remainingIndices{};
  std::ranges::set_difference(std::ranges::iota_view(size_t(0), size_t(rowsFull)), indices, remainingIndices.begin());

  return (E(remainingIndices)).eval();
}

/**
 * \brief Converts a 3x3 matrix to Voigt notation, possibly reducing it based on material properties.
 *  \ingroup utils
 * \tparam ST Scalar type of the matrix.
 * \tparam MaterialImpl Type of the material implementation.
 * \param E Input 3x3 matrix.
 * \param material Reference to the material implementation.
 * \param isStrain Flag indicating if the matrix represents strain (default is true).
 * \return Resulting matrix in Voigt notation.
 *
 * This function converts a 3x3 matrix to its Voigt notation. If the material is not reduced, the full Voigt notation
 * is returned. Otherwise, the specified columns (based on material properties, such as VanishingStress) are removed,
 * and the reduced Voigt notation is returned.
 */
template <typename ST, typename MaterialImpl>
auto toVoigtAndMaybeReduce(const Eigen::Matrix<ST, 3, 3>& E, [[maybe_unused]] const MaterialImpl& material,
                           bool isStrain = true) {
  if constexpr (!MaterialImpl::isReduced)
    return toVoigt(E, isStrain);
  else {
    auto ev = toVoigt(E, isStrain);
    static_assert(decltype(ev)::RowsAtCompileTime == 6 and decltype(ev)::ColsAtCompileTime == 1);

    auto evRed = removeCol(ev, MaterialImpl::fixedVoigtIndices);
    static_assert(decltype(evRed)::RowsAtCompileTime == 6 - MaterialImpl::fixedVoigtIndices.size() and
                  decltype(evRed)::ColsAtCompileTime == 1);
    return evRed;
  }
}

/**
 * \brief Enlarges a matrix if it reduced in the context of material laws, i.e., VanishingStress
 * If the material is not reduced the untouched matrix is returned and rendering the function as a NoOp.
 *  \ingroup utils
 * \tparam M Type of the material.
 * \tparam Derived Type of the input matrix.
 * \param E Input matrix.
 * \return auto Resulting matrix based on material properties.
 * \details
 * This function takes an input matrix and, based on the material properties, either returns the original matrix
 * (if it is not reduced) or enlarges the matrix by filling in the specified columns with zeros (if it is reduced).
 */
template <typename M, typename Derived>
decltype(auto) enlargeIfReduced(const Eigen::MatrixBase<Derived>& E) {
  if constexpr (!Concepts::EigenVector6<Derived> and Concepts::EigenVector<Derived>) {
    static_assert(M::isReduced, "You should only end up here, if your material is reduced");

    auto freeindices = M::MaterialImpl::freeVoigtIndices;
    auto p_E         = Eigen::Vector<typename M::MaterialImpl::ScalarType, 6>::Zero().eval();
    for (int ri = 0; auto i : freeindices)
      p_E(i) = E(ri++);
    return p_E;

  } else if constexpr (Concepts::EigenMatrix66<Derived> or Concepts::EigenMatrix33<Derived> or
                       Concepts::EigenVector6<Derived>) {
    return E.derived();
  } else {
    static_assert(M::isReduced, "You should only end up here, if your material is reduced");

    auto freeindices = M::MaterialImpl::freeVoigtIndices;
    auto p_E         = Eigen::Matrix<typename M::MaterialImpl::ScalarType, 6, 6>::Zero().eval();
    for (int ri = 0; auto i : freeindices) {
      for (int rj = 0; auto j : freeindices)
        p_E(i, j) = E(ri, rj++);
      ++ri;
    }
    return p_E;
  }
}

} // namespace Ikarus
