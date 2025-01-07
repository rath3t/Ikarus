// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file eigensparseaddon.hh
 * \brief Enhance the eigen sparse matrix types with given given functions.
 */

#pragma once

/**
 * \brief Get the linear index corresponding to the given row and column indices.
 * \param row The row index.
 * \param col The column index.
 * \return The linear index.
 * @note The function assumes that `IsRowMajor`, `rows()`, `cols()`, `m_outerIndex`, `m_innerNonZeros`,
 * and `m_data` are accessible within the scope.
 */
inline Index getLinearIndex(Index row, Index col) const {
  eigen_assert(row >= 0 && row < rows() && col >= 0 && col < cols());

  const Index outer = IsRowMajor ? row : col;
  const Index inner = IsRowMajor ? col : row;

  Index start = m_outerIndex[outer];
  Index end   = m_innerNonZeros ? m_outerIndex[outer] + m_innerNonZeros[outer] : m_outerIndex[outer + 1];
  eigen_assert(end >= start && "you probably called coeffRef on a non finalized matrix");

  Index p = m_data.searchLowerIndex(start, end - 1, StorageIndex(inner));
  return p;
}
