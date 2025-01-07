// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file math.hh
 * \brief Implementation of math related algorithms
 */

#pragma once
#include <concepts>

namespace Ikarus {
namespace Impl {
  /**
   * \brief Helper function for compile-time square root calculation.
   *
   * \tparam T The type for which square root is calculated.
   * \param x The value for which square root is calculated.
   * \param lo Lower bound for the search interval.
   * \param hi Upper bound for the search interval.
   * \return constexpr T The calculated square root.
   */
  template <typename T>
  constexpr T sqrt_helper(T x, T lo, T hi) {
    if (lo == hi)
      return lo;

    const T mid = (lo + hi + 1) / 2;

    if (x / mid < mid)
      return sqrt_helper<T>(x, lo, mid - 1);
    else
      return sqrt_helper(x, mid, hi);
  }

} // namespace Impl

/**
 * \brief Compile-time square root for integer types.
 *
 * \tparam T The integral type for which square root is calculated.
 * \param x The value for which square root is calculated.
 * \return constexpr T The calculated square root.
 */
template <typename T>
requires std::integral<T>
constexpr T ct_sqrt(T x) {
  return Impl::sqrt_helper<T>(x, 0, x / 2 + 1);
}
} // namespace Ikarus
