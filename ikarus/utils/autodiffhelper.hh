// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file autodiffhelper.hh
 * \brief Helper for the autodiff library
 */

#pragma once
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace Ikarus::utils {

/**
 * \brief Computes the Hessian matrix for each parameter of a given function.
 * \ingroup  utils
 * The Hessian matrix represents the second-order partial derivatives of the function with respect to the specified
 * variables.
 *
 * \tparam Fun The type of the function to be differentiated.
 * \tparam Vars The types of the variables with respect to which the Hessian is computed.
 * \tparam Args The types of the arguments passed to the function.
 * \tparam U The type representing the result of the function evaluation.
 * \tparam G The type representing the gradient of the function.
 * \tparam H The type representing the Hessian matrix.
 * \param f The function to be differentiated.
 * \param wrt The variables with respect to which the Hessian is computed.
 * \param at The values at which the Hessian is evaluated.
 * \param u The result of the function evaluation.
 * \param g The gradient of the function.
 * \param h The Hessian matrix (output).
 */
template <typename Fun, typename... Vars, typename... Args, typename U, typename G, typename H>
void hessianN(const Fun& f, const autodiff::Wrt<Vars...>& wrt, const autodiff::At<Args...>& at, U& u,
              std::array<G, U::RowsAtCompileTime>& g, std::array<H, U::RowsAtCompileTime>& h) {
  static_assert(sizeof...(Vars) >= 1);
  static_assert(sizeof...(Args) >= 1);

  auto fEntry = [&](auto& I) { return [&](const auto&) { return std::apply(f, at.args)[I]; }; };
  for (int i = 0; i < U::RowsAtCompileTime; ++i)
    hessian(fEntry(i), wrt, at, u[i], g[i], h[i]);
}
} // namespace Ikarus::utils
