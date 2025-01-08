// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file tensorutils.hh
 * \brief Helper for the Eigen::Tensor types
 */

#pragma once

#include <numeric>
#include <ranges>
#include <unsupported/Eigen/CXX11/Tensor>

#include <dune/common/promotiontraits.hh>

#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/math.hh>
namespace Ikarus {

/**
 * \brief View an Eigen matrix as an Eigen Tensor with specified dimensions.
 * \ingroup tensor
 * \tparam Derived Type of the input Eigen matrix.
 * \tparam T Type of the elements in the matrix.
 * \tparam rank Rank of the resulting Tensor.
 * \param matrix Input Eigen matrix to be cast.
 * \param dims Dimensions of the resulting Tensor.
 * \return Eigen::Tensor<typename Derived::Scalar, rank> The casted Eigen Tensor.
 */
template <typename Derived, typename T, auto rank>
Eigen::Tensor<typename Derived::Scalar, rank> tensorView(const Eigen::EigenBase<Derived>& matrix,
                                                         const std::array<T, rank>& dims) {
  return Eigen::TensorMap<const Eigen::TensorFixedSize<
      const typename Derived::Scalar, Eigen::Sizes<Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>>>(
      matrix.derived().eval().data(), dims);
}

/**
 * \brief Computes the dyadic product of two Eigen tensors.
 * \details The components of the result read  \f[ \CI_{ijkl} = A_{ij}B_{kl}. \f]
 * \ingroup tensor
 * \param A_ij First tensor.
 * \param B_kl Second tensor.
 * \return  Resulting tensor after the dyadic product.
 */
auto dyadic(const auto& A_ij, const auto& B_kl) {
  Eigen::array<Eigen::IndexPair<long>, 0> empty_index_list = {};
  return A_ij.contract(B_kl, empty_index_list).eval();
}

/**
 * \brief Computes the dyadic product of two first order Tensors (here: Eigen::Vector).
 * \details The components of the result read  \f[ A_{ij} = a_{i}b_{j}. \f]
 * \ingroup tensor
 * \param a_i First tensor.
 * \param b_j Second tensor.
 * \return  Resulting tensor after the dyadic product
 */
template <typename ST_, int size, bool asTensor = true>
auto dyadic(const Eigen::Vector<ST_, size>& a, const Eigen::Vector<ST_, size>& b) {
  if constexpr (asTensor)
    return tensorView((a * b.transpose()).eval(), std::array<Eigen::Index, 2>({size, size}));
  else
    return (a * b.transpose()).eval();
}

/**
 * \brief Generates a symmetric identity fourth-order tensor.
 * \ingroup tensor
 * \tparam ScalarType Type of the elements in the tensor.
 * \tparam dim Dimension of the tensor.
 * \return  Symmetric identity fourth-order tensor.
 */
template <typename ScalarType = double, int dim = 3>
auto symmetricIdentityFourthOrder() {
  Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>> idTensor;
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      for (int k = 0; k < dim; ++k)
        for (int l = 0; l < dim; ++l)
          idTensor(i, j, k, l) = 0.5 * ((i == k) * (j == l) + (i == l) * (j == k));
  return idTensor;
}

/**
 * \brief Generates a symmetric fourth-order tensor based on two second-order input tensors.
 * \ingroup tensor
 * \details The components of the result read  \f[ \CI_{ijkl} = \frac{1}{2} \left(A_{ik}B_{jl}+A_{il}B_{jk} \right) .
 * \f]
 * \tparam ScalarType Type of the elements in the tensors.
 * \tparam dim Dimension of the tensors.
 * \param A First tensor.
 * \param B Second tensor.
 * \return  Symmetric fourth-order tensor.
 */
template <typename ScalarType = double, int dim = 3>
auto symmetricFourthOrder(const auto& A, const auto& B) {
  Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>> res;
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      for (int k = 0; k < dim; ++k)
        for (int l = 0; l < dim; ++l)
          res(i, j, k, l) = 0.5 * (A(i, k) * B(j, l) + A(i, l) * B(j, k));
  return res;
}

/**
 * \brief Generates an identity fourth-order tensor.
 *  \ingroup tensor
 * \details The components of the result read  \f[ \CI_{ijkl} = \de_{ij}\de_{kl}. \f]
 * \tparam ScalarType Type of the elements in the tensor.
 * \tparam dim Dimension of the tensor.
 * \return  Identity fourth-order tensor.
 */
template <typename ScalarType = double, int dim = 3>
auto identityFourthOrder() {
  Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>> idTensor;
  idTensor.setZero();
  for (int i = 0; i < dim; ++i)
    for (int k = 0; k < dim; ++k)
      idTensor(i, i, k, k) = 1.0;
  return idTensor;
}

/**
 * \brief Computes the IKJL product of two matrices.
 *  \ingroup tensor
 * \details The components of the result read  \f[ \CI_{ijkl} = A_{ik}B_{jl}, \f] which simply swaps the inner slots
 * `j`and `k`
 * \tparam AType Type of the first matrix.
 * \tparam BType Type of the second matrix.
 * \param A First matrix.
 * \param B Second matrix.
 * \return  Resulting tensor of the IKJL product.
 */
template <typename AType, typename BType>
auto fourthOrderIKJL(const Eigen::MatrixBase<AType>& A, const Eigen::MatrixBase<BType>& B) {
  static_assert(AType::RowsAtCompileTime == BType::RowsAtCompileTime);
  static_assert(AType::ColsAtCompileTime == BType::ColsAtCompileTime);
  using ScalarResultType = typename Dune::PromotionTraits<typename AType::Scalar, typename BType::Scalar>::PromotedType;
  constexpr int dim      = AType::RowsAtCompileTime;
  Eigen::TensorFixedSize<ScalarResultType, Eigen::Sizes<dim, dim, dim, dim>> res;
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      for (int k = 0; k < dim; ++k)
        for (int l = 0; l < dim; ++l)
          res(i, j, k, l) = A(i, k) * B(j, l);
  return res;
}

/**
 * \brief Creates a symmetric fourth-order tensor in the two specified slots of the input tensor.
 *  \ingroup tensor
 * \tparam ScalarType Type of the elements in the tensor.
 * \param t Input tensor.
 * \param slots Indices of the slots to be swapped.
 * \return  Fourth-order Tensor which is symmetric in the given slots.
 */
template <typename ScalarType, long int dim>
auto symTwoSlots(const Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>>& t,
                 const std::array<size_t, 2>& slots) {
  std::array<size_t, 4> shuffleSlot;
  std::iota(shuffleSlot.begin(), shuffleSlot.end(), 0);    // create 0,1,2,3 array
  std::swap(shuffleSlot[slots[0]], shuffleSlot[slots[1]]); // swap  the given slots
  return (0.5 * (t + t.shuffle(shuffleSlot))).eval();
}

/**
 * \brief Converts 2D indices to Voigt notation index.
 *  \ingroup tensor
 * \param i Row index.
 * \param j Column index.
 * \return Eigen::Index Voigt notation index.
 *
 * This function converts 2D indices (i, j) to a Voigt notation index.
 * The Voigt notation is used to represent the six unique components of a symmetric 3x3 matrix
 * in a one-dimensional array.
 *
 * If the input indices are not within the valid range (0, 1, 2), an assertion failure is triggered.
 */
constexpr Eigen::Index toVoigt(Eigen::Index i, Eigen::Index j) noexcept {
  if (i == j) // _00 -> 0, _11 -> 1,  _22 -> 2
    return i;
  if ((i == 1 and j == 2) or (i == 2 and j == 1)) // _12 and _21 --> 3
    return 3;
  if ((i == 0 and j == 2) or (i == 2 and j == 0)) // _02 and _20 --> 4
    return 4;
  if ((i == 0 and j == 1) or (i == 1 and j == 0)) // _01 and _10 --> 5
    return 5;
  assert(i < 3 and j < 3 && "For Voigt notation the indices need to be 0,1 or 2.");
  __builtin_unreachable();
}

/**
 * \brief Converts a fourth-order tensor of fixed size 3x3x3x3 to a Voigt notation matrix of size 6x6.
 *  \ingroup tensor
 * \tparam ScalarType Data type of the tensor elements.
 * \param ft Fourth-order tensor .
 * \return Voigt notation matrix.
 *
 * This function converts a fourth-order tensor to a Voigt notation matrix, which is a symmetric 6x6 matrix
 * containing the unique components of the input tensor. The mapping from the tensor indices to the Voigt notation
 * indices is performed by the toVoigt function.
 *
 * \remarks The current implementation
 * does not take advantage of this symmetry.
 */
template <typename ScalarType = double>
Eigen::Matrix<ScalarType, 6, 6> toVoigt(const Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>>& ft) {
  Eigen::Matrix<ScalarType, 6, 6> mat;
  for (Eigen::Index i = 0; i < 3; ++i)
    for (Eigen::Index j = 0; j < 3; ++j)
      for (Eigen::Index k = 0; k < 3; ++k)
        for (Eigen::Index l = 0; l < 3; ++l)
          mat(toVoigt(i, j), toVoigt(k, l)) = ft(i, j, k, l);
  return mat;
}

/**
 * \brief Converts a square 2x2 or 3x3 matrix to a Voigt notation vector.
 *  \ingroup tensor
 * \tparam ST Data type of the matrix elements.
 * \tparam size Number of rows and columns of the square matrix.
 * \tparam Options Eigen matrix options.
 * \param E Input matrix of size (size x size).
 * \param isStrain Flag indicating whether the conversion is for strain (true) or not (false).
 * \return  Vector with components in Voigt notation vector.
 *
 * This function converts a square matrix to a Voigt notation vector, which contains the unique components of
 * the input matrix. The mapping from the matrix indices to the Voigt notation indices is performed by the toVoigt
 * function.
 *
 * The optional isStrain parameter allows the user to specify whether the conversion is intended for strain
 * calculations. If isStrain is true, the off-diagonal components are multiplied by 2, providing the correct Voigt
 * notation for symmetric strain tensors.
 */

template <typename ST, int size, int Options, int maxSize>
requires((size > 0 and size <= 3) or (maxSize > 0 and maxSize <= 3 and size == Eigen::Dynamic))
auto toVoigt(const Eigen::Matrix<ST, size, size, Options, maxSize, maxSize>& E, bool isStrain = true) {
  constexpr bool isFixedSized   = (size != Eigen::Dynamic);
  const ST possibleStrainFactor = isStrain ? 2.0 : 1.0;

  const size_t inputSize = isFixedSized ? size : E.rows();
  auto EVoigt            = [&]() {
    if constexpr (isFixedSized) {
      Eigen::Vector<ST, (size * (size + 1)) / 2> EVoigt;
      EVoigt.template head<size>() = E.diagonal();
      return EVoigt;
    } else {
      Eigen::Matrix<ST, Eigen::Dynamic, 1, Options, 6, 1> EVoigt;
      EVoigt.resize((inputSize * (inputSize + 1)) / 2);
      EVoigt.template head(inputSize) = E.diagonal();
      return EVoigt;
    }
  }();

  if (inputSize == 2)
    EVoigt(2) = E(0, 1) * possibleStrainFactor;
  else if (inputSize == 3) {
    EVoigt(inputSize)     = E(1, 2) * possibleStrainFactor;
    EVoigt(inputSize + 1) = E(0, 2) * possibleStrainFactor;
    EVoigt(inputSize + 2) = E(0, 1) * possibleStrainFactor;
  }
  return EVoigt;
}

/**
 * \brief Converts a vector given in Voigt notation to a matrix.
 *  \ingroup tensor
 * \tparam ST Scalar type of the vector elements.
 * \tparam size Size of the Voigt notation vector.
 * \param EVoigt Voigt notation vector.
 * \param isStrain Flag indicating whether the vector represents a strain (default is true).
 * \return Matrix corresponding to the vector in Voigt notation.
 *  \details
 * This function converts a vector given in Voigt notation to the corresponding matrix. The conversion depends on the
 * size The parameter `isStrain` is used to determine the conversion factor for off-diagonal components, which need to
 * be divided by 2 in the matrix representation if the quantity is a strain tensor.
 *
 * The function requires that the size of the Voigt notation vector is valid (1, 3, or 6).
 */
template <typename ST, int size, int Options, int maxSize>
requires((size == 1 or size == 3 or size == 6) or
         ((maxSize == 1 or maxSize == 3 or maxSize == 6) and size == Eigen::Dynamic))
auto fromVoigt(const Eigen::Matrix<ST, size, 1, Options, maxSize, 1>& EVoigt, bool isStrain = true) {
  constexpr bool isFixedSized   = (size != Eigen::Dynamic);
  const ST possibleStrainFactor = isStrain ? 0.5 : 1.0;

  const size_t inputSize = isFixedSized ? size : EVoigt.size();
  const size_t matrixSize =
      isFixedSized ? (-1 + ct_sqrt(1 + 8 * size)) / 2 : (-1 + static_cast<int>(std::sqrt(1 + 8 * inputSize))) / 2;

  auto E = [&]() {
    if constexpr (isFixedSized) {
      Eigen::Matrix<ST, matrixSize, matrixSize> E;
      E.diagonal() = EVoigt.template head<matrixSize>();
      return E;
    } else {
      Eigen::Matrix<ST, Eigen::Dynamic, Eigen::Dynamic, Options, 3, 3> E;
      E.resize(matrixSize, matrixSize);
      E.diagonal() = EVoigt.template head(matrixSize);
      return E;
    }
  }();

  if (matrixSize == 2) {
    E(0, 1) = E(1, 0) = EVoigt(2) * possibleStrainFactor;
  } else if (matrixSize == 3) {
    E(2, 1) = E(1, 2) = EVoigt(matrixSize) * possibleStrainFactor;
    E(2, 0) = E(0, 2) = EVoigt(matrixSize + 1) * possibleStrainFactor;
    E(1, 0) = E(0, 1) = EVoigt(matrixSize + 2) * possibleStrainFactor;
  }

  return E;
}

/**
 * \brief Converts a Voigt notation index to matrix indices.
 *  \ingroup tensor
 * \param i Voigt notation index.
 * \return Matrix indices corresponding to the Voigt notation index.
 *  \details
 * This function converts a Voigt notation index to the corresponding matrix indices. The mapping is based on the
 * assumption that the Voigt notation indices 0, 1, and 2 represent the diagonal components `00`, `11`, and `22`,
 * respectively. The remaining Voigt notation indices (3, 4, and 5) correspond to the off-diagonal components
 * (`12` and `21`, `02` and `20`, `01` and `10`).
 *
 * The function asserts that the input index is within the valid range for Voigt notation (0 to 5).
 */
constexpr std::array<size_t, 2> fromVoigt(size_t i) {
  if (i < 3) // _00 -> 0, _11 -> 1,  _22 -> 2
    return {i, i};
  else if (i == 3)
    return {1, 2};
  else if (i == 4)
    return {0, 2};
  else if (i == 5)
    return {0, 1};
  else {
    assert(i < 6 && "For Voigt notation the indices need to be 0 and 5.");
    __builtin_unreachable();
  }
}

/**
 * \brief Converts a matrix in Voigt notation to a Fourth-order tensor.
 *  \ingroup tensor
 * \tparam ScalarType Scalar type of the matrix elements.
 * \param CVoigt Voigt notation matrix.
 * \return Fourth-order tensor corresponding to the matrix in Voigt notation.
 *  \details
 * This function converts a Voigt notation matrix to the corresponding 4th-order tensor. The function uses the
 * `fromVoigt` function to map matrix indices to tensor indices. The resulting tensor is symmetric due to symmetry
 * considerations.
 */
template <typename ScalarType>
auto fromVoigt(const Eigen::Matrix<ScalarType, 6, 6>& CVoigt) {
  Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> C;
  // size_t iR=0,jR=0;
  for (size_t i = 0; i < 6; ++i) {
    for (size_t j = 0; j < 6; ++j) {
      auto firstIndices                                                       = fromVoigt(i);
      auto secondIndices                                                      = fromVoigt(j);
      C(firstIndices[0], firstIndices[1], secondIndices[0], secondIndices[1]) = CVoigt(i, j);
    }
  }
  return C;
}

/**
 * \brief Calculates the 2D transformation matrix.
 *
 * \details This function computes the transformation matrix needed to transform second-order tensors
 * represented in Voigt notation from local to global coordinate system for 2D elements.
 *
 * \tparam Geometry The geometry type.
 * \param geometry Reference to the geometry object.
 * \return The transformation matrix for 2D elements.
 */
template <typename Geometry>
Eigen::Matrix3d calcTransformationMatrix2D(const Geometry& geometry) {
  const auto& referenceElement = Dune::ReferenceElements<double, 2>::general(geometry.type());
  const auto quadPos0          = referenceElement.position(0, 0);

  const auto jacobianinvT0 = toEigen(geometry.jacobianInverseTransposed(quadPos0));
  const auto detJ0         = geometry.integrationElement(quadPos0);

  auto jaco = (jacobianinvT0).inverse().eval();
  auto J11  = jaco(0, 0);
  auto J12  = jaco(0, 1);
  auto J21  = jaco(1, 0);
  auto J22  = jaco(1, 1);

  Eigen::Matrix3d T0{
      {      J11 * J11,       J12 * J12,             J11 * J12},
      {      J21 * J21,       J22 * J22,             J21 * J22},
      {2.0 * J11 * J21, 2.0 * J12 * J22, J21 * J12 + J11 * J22}
  };

  return T0.inverse() * detJ0;
}

/**
 * \brief Calculates the 3D transformation matrix.
 *
 * \details This function computes the transformation matrix needed to transform second-order tensors
 * represented in Voigt notation from local to global coordinate system for 3D elements.
 *
 * \tparam Geometry The geometry type.
 * \param geometry Reference to the geometry object.
 * \return The transformation matrix for 3D elements.
 */
template <typename Geometry>
Eigen::Matrix<double, 6, 6> calcTransformationMatrix3D(const Geometry& geometry) {
  const auto& referenceElement = Dune::ReferenceElements<double, 3>::general(geometry.type());
  const auto quadPos0          = referenceElement.position(0, 0);

  const auto jacobianinvT0 = toEigen(geometry.jacobianInverseTransposed(quadPos0));
  const auto detJ0         = geometry.integrationElement(quadPos0);

  auto jaco = (jacobianinvT0).inverse().eval();
  auto J11  = jaco(0, 0);
  auto J12  = jaco(0, 1);
  auto J13  = jaco(0, 2);
  auto J21  = jaco(1, 0);
  auto J22  = jaco(1, 1);
  auto J23  = jaco(1, 2);
  auto J31  = jaco(2, 0);
  auto J32  = jaco(2, 1);
  auto J33  = jaco(2, 2);

  // clang-format off
  Eigen::Matrix<double, 6, 6> T0  {
    {J11 * J11,       J12 * J12,       J13 * J13,             J11 * J12,             J11 * J13,             J12 * J13},
    {J21 * J21,       J22 * J22,       J23 * J23,             J21 * J22,             J21 * J23,             J22 * J23},
    {J31 * J31,       J32 * J32,       J33 * J33,             J31 * J32,             J31 * J33,             J32 * J33},
    {2.0 * J11 * J21, 2.0 * J12 * J22, 2.0 * J13 * J23, J11 * J22 + J21 * J12, J11 * J23 + J21 * J13, J12 * J23 + J22 * J13},
    {2.0 * J11 * J31, 2.0 * J12 * J32, 2.0 * J13 * J33, J11 * J32 + J31 * J12, J11 * J33 + J31 * J13, J12 * J33 + J32 * J13},
    {2.0 * J31 * J21, 2.0 * J32 * J22, 2.0 * J33 * J23, J31 * J22 + J21 * J32, J31 * J23 + J21 * J33, J32 * J23 + J22 * J33}
  };
  // clang-format on

  return T0.inverse() * detJ0;
}

} // namespace Ikarus
