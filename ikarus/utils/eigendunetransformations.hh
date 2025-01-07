// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file eigendunetransformations.hh
 * \brief Helper for transform between Dune linear algebra types and Eigen
 */

#pragma once

#include <assert.h>

#include <dune/common/diagonalmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/scaledidmatrix.hh>

#include <Eigen/Core>

namespace Ikarus {

/**
 * \brief Create Eigen::Vector to Dune::FieldVector.
 * \tparam ST The scalar type of the vectors.
 * \tparam size The size of the vectors.
 * \param vec The Eigen::Vector to be converted.
 * \return Dune::FieldVector<ST, size> representing the converted vector.
 */
template <typename ST, int size>
Dune::FieldVector<ST, size> toDune(const Eigen::Vector<ST, size>& vec) {
  Dune::FieldVector<ST, size> fieldVector;
  for (int i = 0; i < size; ++i)
    fieldVector[i] = vec[i];
  return fieldVector;
}

/**
 * \brief Convert Eigen::Matrix to Dune::FieldVector.
 * \details The matrix has fixed row size. The matrix needs
 * to have a single column.
 * \tparam ST The scalar type of the vectors.
 * \tparam rows The number of rows in the matrix.
 * \param vec The Eigen::Matrix to be converted.
 * \return Dune::FieldVector<ST, rows> representing the converted matrix.
 */
template <typename ST, int rows>
Dune::FieldVector<ST, rows> toDune(const Eigen::Matrix<ST, rows, 0>& vec) {
  assert(vec.cols() == 1 && "The passed matrix needs to have a single column.");
  Dune::FieldVector<ST, rows> fieldVector{0.0};

  for (int i = 0; i < vec.rows(); ++i)
    fieldVector[i] = vec(i, 0);
  return fieldVector;
}

/**
 * \brief Convert Eigen::Matrix to Dune::FieldMatrix.
 * \details The matrix has fixed rows and column size
 * \tparam ST The scalar type of the matrix.
 * \tparam rows The number of rows in the matrix.
 * \tparam cols The number of columns in the matrix.
 * \param mat The Eigen::Matrix to be converted.
 * \return Dune::FieldMatrix<ST, rows, cols> representing the converted matrix.
 */
template <typename ST, int rows, int cols>
Dune::FieldMatrix<ST, rows, cols> toDune(const Eigen::Matrix<ST, rows, cols>& mat) {
  Dune::FieldMatrix<ST, rows, cols> fieldMatrix{0.0};

  for (int i = 0; i < mat.rows(); ++i)
    for (int j = 0; j < mat.cols(); ++j)
      fieldMatrix[i][j] = mat(i, j);
  return fieldMatrix;
}

/**
 * \brief View a Dune::FieldVector as an Eigen::Vector using Map, no copies take place.
 * \tparam ST The scalar type of the vector.
 * \tparam size The size of the vector.
 * \param vec The Dune::FieldVector to be viewed as Eigen::Vector.
 * \return Eigen::Map<const Eigen::Vector<ST, size>> representing the viewed vector.
 */
template <typename ST, int size>
Eigen::Map<const Eigen::Vector<ST, size>> toEigenMap(const Dune::FieldVector<ST, size>& vec) {
  return {vec.data(), size};
}

/**
 * \brief View a constant Dune::FieldVector as a constant Eigen::Vector, no copies take place.
 * \tparam ST The scalar type of the vector.
 * \tparam size The size of the vector.
 * \param vec The Dune::FieldVector to be viewed as Eigen::Vector.
 * \return Eigen::Map<Eigen::Vector<ST, size>> representing the viewed vector.
 */
template <typename ST, int size>
Eigen::Map<Eigen::Vector<ST, size>> toEigenMap(Dune::FieldVector<ST, size>& vec) {
  return {vec.data(), size};
}

} // namespace Ikarus
