//
// Created by Alex on 10.05.2021.
//

#pragma once

#include <dune/common/diagonalmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <Eigen/Core>

namespace Ikarus {

  template <typename ScalarType, int size>
  Dune::FieldVector<ScalarType, size> toFieldVector(const Eigen::Vector<ScalarType, size>& vec) {
    Dune::FieldVector<ScalarType, size> fieldvec;
    for (int i = 0; i < size; ++i)
      fieldvec[i] = vec[i];
    return fieldvec;
  }

  /** \brief Views a dune fieldvector as an Eigen::Vector as Map, no copies take place! */
  template <typename ScalarType, int size>
  Eigen::Map<const Eigen::Vector<ScalarType, size>> toEigenVectorMap(const Dune::FieldVector<ScalarType, size>& vec) {
    return {vec.data(), size};
  }

  template <typename ScalarType, int size>
  Eigen::Vector<ScalarType, size> toEigenVector(const Dune::FieldVector<ScalarType, size>& vec) {
    Eigen::Vector<ScalarType, size> eigenVector;
    for (int i = 0; i < size; ++i)
      eigenVector(i) = vec[i];
    return eigenVector;
  }

  /** \brief Views a const dune fieldvector as a const Eigen::Vector, no copies take place! */
  template <typename ScalarType, int size>
  Eigen::Map<Eigen::Vector<ScalarType, size>> toEigenVector(Dune::FieldVector<ScalarType, size>& vec) {
    return {vec.data(), size};
  }

  template <typename ScalarType, int size1, int size2>
  Eigen::Matrix<ScalarType, size1, size2> toEigenMatrix(const Dune::FieldMatrix<ScalarType, size1, size2>& mat) {
    Eigen::Matrix<ScalarType, size1, size2> eigenmatrix;
    for (int i = 0; i < size1; ++i)
      for (int j = 0; j < size2; ++j)
        eigenmatrix(i, j) = mat[i][j];
    return eigenmatrix;
  }

  template <typename ScalarType, int size1>
  Eigen::Matrix<ScalarType, size1, size1> toEigenMatrix(const Dune::DiagonalMatrix<ScalarType, size1>& mat) {
    Eigen::Matrix<ScalarType, size1, size1> eigenmatrix;
    eigenmatrix.setZero();
    for (int i = 0; i < size1; ++i)
      eigenmatrix(i, i) = mat[i][i];
    return eigenmatrix;
  }

}  // namespace Ikarus