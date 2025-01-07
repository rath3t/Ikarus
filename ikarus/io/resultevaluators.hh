// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file resultevaluators.hh
 * \brief Ikarus Result Evaluators for special stress quantities
 * \ingroup resultevaluators
 *
 */

#pragma once

#include <dune/common/math.hh>

#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::ResultEvaluators {

/**
 * \brief Struct for calculating von Mises stress
 * \ingroup resultevaluators
 * \details The VonMises struct provides a function call operator to calculate von Mises stress.
 * In 2D, this assumes a plane stress state
 */
struct VonMises
{
  /**
   * \brief Calculate the result quantity (von Mises stress)
   * \param resultArray EigenMatrix containing the stress state
   * \param comp component of result (not used here)
   * \tparam R Type of the matrix
   * \return von Mises stress
   */
  template <typename R>
  double operator()(const R& resultArray, [[maybe_unused]] const int comp) const {
    const auto s_x = resultArray(0, 0);
    const auto s_y = resultArray(1, 0);
    if constexpr (R::CompileTimeTraits::RowsAtCompileTime == 3) {
      const auto s_xy = resultArray(2, 0);
      return std::sqrt(Dune::power(s_x, 2) + Dune::power(s_y, 2) - s_x * s_y + 3 * Dune::power(s_xy, 2));
    } else {
      const auto s_z  = resultArray(2, 0);
      const auto s_yz = resultArray(3, 0);
      const auto s_xz = resultArray(4, 0);
      const auto s_xy = resultArray(5, 0);

      return std::sqrt(Dune::power(s_x, 2) + Dune::power(s_y, 2) + Dune::power(s_z, 2) - s_x * s_y - s_x * s_z -
                       s_y * s_z + 3 * (Dune::power(s_xy, 2) + Dune::power(s_xz, 2) + Dune::power(s_yz, 2)));
    }
  }

  /**
   * \brief Get the name of the result type (VonMises)
   * \return String representing the name
   */
  static std::string name() { return "VonMises"; }

  /**
   * \brief Get the number of components in the result (always 1 for VonMises)
   * \return Number of components
   */
  static int ncomps() { return 1; }
};

/**
 * \brief Struct for calculating principal stresses
 * \ingroup resultevaluators
 * \details The PrincipalStress struct provides a function call operator to calculate principal stresses.
 * The components are ordered in a descending manner ($\sigma_1 > \sigma_2$)
 * \tparam dim dimension of stress state
 */
template <int dim>
requires(dim == 2 or dim == 3)
struct PrincipalStress
{
  /**
   * \brief Calculate the result quantity (principal stress)
   * \param resultArray EigenMatrix containing the stress state
   * \param comp component of result
   * \return principal stress
   */
  double operator()(const auto& resultArray, const int comp) const {
    auto mat = fromVoigt(resultArray, false);
    Eigen::SelfAdjointEigenSolver<decltype(mat)> eigensolver(mat, Eigen::EigenvaluesOnly);
    return eigensolver.eigenvalues()[dim - 1 - comp];
  }

  /**
   * \brief Get the name of the result type (PrincipalStress)
   * \return String representing the name
   */
  static std::string name() { return "PrincipalStress"; }

  /**
   * \brief Get the number of components in the result
   * \return Number of components
   */
  static int ncomps() { return dim; }
};

} // namespace Ikarus::ResultEvaluators
