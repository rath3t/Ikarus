// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ranges>

#include <Eigen/Dense>

#include <ikarus/utils/concepts.hh>

namespace Ikarus::Impl {

/**
 * \brief Represents a pair of stress matrix indices (row and column).
 */
struct MatrixIndexPair
{
  Eigen::Index row; ///< Row index.
  Eigen::Index col; ///< Column index.
};

/**
 * \brief Helper function to create an array of free Voigt indices.
 * \tparam size The size of the fixed pairs array.
 * \param fixed An array of MatrixIndexPair representing fixed indices.
 * \return std::array<size_t, 6 - size> The array of free Voigt indices.
 */
template <size_t size>
consteval auto createfreeVoigtIndices(const std::array<MatrixIndexPair, size>& fixed) {
  std::array<size_t, 6 - size> res{};
  std::array<size_t, size> voigtFixedIndices;
  std::ranges::transform(fixed, voigtFixedIndices.begin(), [](auto pair) { return toVoigt(pair.row, pair.col); });
  std::ranges::sort(voigtFixedIndices);
  std::ranges::set_difference(std::ranges::iota_view(size_t(0), size_t(6)), voigtFixedIndices, res.begin());
  std::ranges::sort(res);
  return res;
}

/**
 * \brief Helper function to create an array of fixed Voigt indices.
 * \tparam size The size of the fixed pairs array.
 * \param fixed An array of MatrixIndexPair representing fixed indices.
 * \return std::array<size_t, size> The array of fixed Voigt indices.
 */
template <size_t size>
consteval auto createFixedVoigtIndices(const std::array<MatrixIndexPair, size>& fixed) {
  std::array<size_t, size> fixedIndices;
  std::ranges::transform(fixed, fixedIndices.begin(), [](auto pair) { return toVoigt(pair.row, pair.col); });
  std::ranges::sort(fixedIndices);
  return fixedIndices;
}

/**
 * \brief Helper function to count the number of diagonal indices in the fixed pairs array.
 * \tparam size The size of the fixed pairs array.
 * \param fixed An array of MatrixIndexPair representing fixed indices.
 * \return constexpr size_t The number of diagonal indices.
 */
template <size_t size>
constexpr size_t countDiagonalIndices(const std::array<MatrixIndexPair, size>& fixed) {
  size_t count = 0;
  for (auto v : fixed) {
    if (v.col == v.row)
      ++count;
  }
  return count;
}

/**
 * \brief Converts the input strain matrix to the appropriate form for stress reduction.
 * \tparam Derived The derived type of the input matrix.
 * \param E The input strain matrix.
 * \return decltype(auto) The converted strain matrix.
 */
template <typename Derived>
decltype(auto) maybeFromVoigt(const Eigen::MatrixBase<Derived>& E) {
  if constexpr (Concepts::EigenVector<Derived>) { // receiving vector means Voigt notation
    return fromVoigt(E.derived(), true);
  } else
    return E.derived();
}

} // namespace Ikarus::Impl