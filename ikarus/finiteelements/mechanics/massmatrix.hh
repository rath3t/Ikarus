// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file massmatrix.hh
 * \brief Standard implementation for trivial mass matrices
 * \ingroup  mechanics
 */

#pragma once

#include <unsupported/Eigen/KroneckerProduct>

namespace Ikarus {
/**
 * \brief Evaluates the Kronecker Product to obtain the entries of a mass Matrix.
 *
 * \tparam dim degrees of freedom per node.
 * \tparam Derived the derived Matrix type.
 * \param intElement factor for numerical integration.
 * \param N the shape functions.
 * \param rho density (or any factor).
 * \param M the Mass matrix that the Kronecker Product should be evaluated into.
 */
template <int dim, typename Derived>
void evaluateKroneckerProduct(const auto& intElement, const auto& N, const auto& rho, Eigen::MatrixBase<Derived>& M) {
  M += Eigen::kroneckerProduct((N * N.transpose()).eval(),
                               Eigen::Matrix<double, dim, dim>::Identity() * intElement * rho)
           .eval();
}
} // namespace Ikarus
