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

#include <ikarus/finiteelements/mechanics/materials.hh>
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
   * \param resultArray EigenVector containing the stress state in Voigt notation
   * \param comp component of result (not used here)
   * \tparam R Type of the matrix
   * \return von Mises stress
   */
  template <typename R>
  double operator()(const R& resultArray, [[maybe_unused]] const int comp) const {
    const auto s_x = resultArray[0];
    const auto s_y = resultArray[1];
    if constexpr (R::CompileTimeTraits::RowsAtCompileTime == 3) {
      const auto s_xy = resultArray[2];
      return std::sqrt(Dune::power(s_x, 2) + Dune::power(s_y, 2) - s_x * s_y + 3 * Dune::power(s_xy, 2));
    } else {
      const auto s_z  = resultArray[2];
      const auto s_yz = resultArray[3];
      const auto s_xz = resultArray[4];
      const auto s_xy = resultArray[5];

      return std::sqrt(Dune::power(s_x, 2) + Dune::power(s_y, 2) + Dune::power(s_z, 2) - s_x * s_y - s_x * s_z -
                       s_y * s_z + 3 * (Dune::power(s_xy, 2) + Dune::power(s_xz, 2) + Dune::power(s_yz, 2)));
    }
  }

  /**
   * \brief Get the name of the result type
   * \return String representing the name
   */
  constexpr static std::string name() { return "VonMises"; }

  /**
   * \brief Get the number of components in the result (always 1 for VonMises)
   * \return Number of components
   */
  constexpr static int ncomps() { return 1; }
};

/**
 * \brief Struct for calculating hydrostatic stress
 * \ingroup resultevaluators
 * \details The HydrostaticStress struct provides a function call operator to calculate hydrostatic stress.
 * In 2D, this assumes a plane stress state
 */
struct HydrostaticStress
{
  /**
   * \brief Calculate the result quantity (hydrostatic stress)
   * \param resultArray EigenVector containing the stress state in Voigt notation
   * \param comp component of result (not used here)
   * \tparam R Type of the matrix
   * \return Hydrostatic stress
   */
  template <typename R>
  double operator()(const R& resultArray, [[maybe_unused]] const int comp) const {
    static constexpr int size = R::CompileTimeTraits::RowsAtCompileTime;
    const auto sigma          = fromVoigt(resultArray, false);
    constexpr static int dim  = std::remove_cvref_t<decltype(sigma)>::CompileTimeTraits::RowsAtCompileTime;
    return 1.0 / dim * sigma.trace();
  }

  /**
   * \brief Get the name of the result type
   * \return String representing the name
   */
  constexpr static std::string name() { return "HydrostaticStress"; }

  /**
   * \brief Get the number of components in the result (always 1 for hydrostatic stress)
   * \return Number of components
   */
  constexpr static int ncomps() { return 1; }
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
   * \param resultArray EigenVector containing the stress state in Voigt notation
   * \param comp component of result
   * \return principal stress
   */
  double operator()(const auto& resultArray, const int comp) const {
    auto mat = fromVoigt(resultArray, false);
    Eigen::SelfAdjointEigenSolver<decltype(mat)> eigensolver(mat, Eigen::EigenvaluesOnly);
    return eigensolver.eigenvalues()[dim - 1 - comp];
  }

  /**
   * \brief Get the name of the result type
   * \return String representing the name
   */
  constexpr static std::string name() { return "PrincipalStress"; }

  /**
   * \brief Get the number of components in the result
   * \return Number of components
   */
  constexpr static int ncomps() { return dim; }
};

/**
 * \brief Struct for calculating stress triaxiality
 * \ingroup resultevaluators
 * \details The Triaxiality struct provides a function call operator to calculate stress triaxiality.
 * In 2D, this assumes a plane stress state
 */
struct Triaxiality
{
  /**
   * \brief Calculate the result quantity (stress triaxiality)
   * \param resultArray EigenVector containing the stress state in Voigt notation
   * \param comp component of result (not used here)
   * \tparam R Type of the matrix
   * \return Triaxiality stress
   */
  template <typename R>
  double operator()(const R& resultArray, [[maybe_unused]] const int comp) const {
    auto sigeq = VonMises{}(resultArray, 0);
    auto sigm  = HydrostaticStress{}(resultArray, 0);
    return sigm / sigeq;
  }
  /**
   * \brief Get the name of the result type
   * \return String representing the name
   */
  constexpr static std::string name() { return "Triaxiality"; }

  /**
   * \brief Get the number of components in the result  (always 1 for stress triaxiality)
   * \return Number of components
   */
  constexpr static int ncomps() { return 1; }
};

/**
 * \brief Struct for calculating the 2d polar stress. The center of the coordinate system is chosen to be the center of
 * the corresponding reference geometry.
  \ingroup resultevaluators
 */
struct PolarStress
{
  PolarStress(const Dune::FieldVector<double, 2>& origin)
      : origin_(origin) {}

  /**
   * \brief Calculate the result quantity (von Mises stress)
   * \param resultArray EigenVector containing the stress state in Voigt notation
   * \param comp component of result
   * \tparam R Type of the matrix
   * \return von Mises stress
   */
  template <typename R>
  double operator()(const R& resultArray, const auto& pos, const auto& fe, const int comp) const {
    static_assert(R::CompileTimeTraits::RowsAtCompileTime == 3, "PolarStress is only valid for 2D.");

    // Offset to center the coordinate system in the reference geometry
    Dune::FieldVector<double, 2> posGlobal = fe.geometry().global(pos) - origin_;
    auto theta                             = std::atan2(posGlobal[1], posGlobal[0]);

    const auto s_x  = resultArray[0];
    const auto s_y  = resultArray[1];
    const auto s_xy = resultArray[2];

    auto hpl   = 0.5 * (s_x + s_y);
    auto hmi   = 0.5 * (s_x - s_y);
    auto cos2t = std::cos(2 * theta);
    auto sin2t = std::sin(2 * theta);

    switch (comp) {
      case 0:
        return hpl + hmi * cos2t + s_xy * sin2t;
      case 1:
        return hpl - hmi * cos2t - s_xy * sin2t;
      case 2:
        return -hmi * sin2t + s_xy * cos2t;
      default:
        DUNE_THROW(Dune::RangeError, "PolarStress: Comp out of range.");
    }
  }

  /**
   * \brief Get the name of the result type
   * \return String representing the name
   */
  constexpr static std::string name() { return "PolarStress"; }

  /**
   * \brief Get the number of components in the result
   * \return Number of components
   */
  constexpr static int ncomps() { return 3; }

private:
  Dune::FieldVector<double, 2> origin_;
};

} // namespace Ikarus::ResultEvaluators
