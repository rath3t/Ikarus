// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file lumpingschemes.hh
 * \brief Implementation of lumping schemes
 */

#pragma once
#include <dune/common/rangeutilities.hh>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Ikarus::Dynamics::LumpingSchemes {

/**
 * \brief Implements the operator() for a row-sum lumping scheme that can be passed to the assemblermanipulator to lump
 * the already assembled global mass matrix.
 */
struct RowSumLumping
{
  /**
   * \brief Implements row-sum lumping for sparse matrices. Warning: the new zero-entries are getting pruned (i.e.
   * deleted)
   *
   * \tparam ScalarType Scalartype of the sparse matrix
   * \param mat the sparse matrix to be modified
   */
  template <typename ScalarType = double>
  void operator()(const auto&, const auto&, auto, auto, Eigen::SparseMatrix<ScalarType>& mat) {
    for (auto i : Dune::range(mat.rows())) {
      auto sum           = mat.row(i).sum();
      mat.coeffRef(i, i) = sum;
    }
    // Deletes all entries expect for main diagonal (entriese are really getting deleted)
    mat.prune([](int i, int j, auto) { return i == j; });
  }

  /**
   * \brief Implements row-sum lumping for dense matrices.
   *
   * \tparam ScalarType Scalartype of the dense matrix
   * \param mat the dense matrix to be modified
   */
  template <typename ScalarType = double>
  void operator()(const auto&, const auto&, auto, auto, Eigen::MatrixX<ScalarType>& mat) {
    for (auto i : Dune::range(mat.rows())) {
      auto sum  = mat.row(i).sum();
      mat(i, i) = sum;
    }
    mat = mat.diagonal().asDiagonal();
  }
};

} // namespace Ikarus::Dynamics::LumpingSchemes
