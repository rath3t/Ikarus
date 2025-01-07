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

#include <ikarus/finiteelements/mechanics/materials/vanishingstrain.hh>
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
   * \brief Get the name of the result type (VonMises)
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
 * \brief Struct for calculating von Mises stress
 * \ingroup resultevaluators
 * \details The VonMises struct provides a function call operator to calculate von Mises stress.
 * In 2D, this assumes a plane stress state
 */
struct HydrostaticStress
{
  /**
   * \brief Calculate the result quantity (von Mises stress)
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
   * \brief Get the number of components in the result (always 1 for HydrostaticStress)
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
   * \brief Get the name of the result type (PrincipalStress)
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
 * \brief Struct for calculating Triaxiality stresses
 * \ingroup resultevaluators
 * \details The VonMises struct provides a function call operator to calculate von Mises stress.
 * In 2D, this assumes a plane stress state
 */
struct Triaxiality
{
  /**
   * \brief Calculate the result quantity (von Mises stress)
   * \param resultArray EigenVector containing the stress state in Voigt notation
   * \param comp component of result (not used here)
   * \tparam R Type of the matrix
   * \return Triaxiality stress
   */
  template <typename R>
  double operator()(const R& resultArray, const int comp) const {
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
   * \brief Get the number of components in the result  (always 1 for Triaxiality)
   * \return Number of components
   */
  constexpr static int ncomps() { return 1; }
};

/**
 * \brief Wrapper to obtain stress results for vanishing materials. It takes a resultevaluator as template argument.
 * For now this works only for plane strain case.
 * If you just want to obtain the 3d stresses without using a resultevaluator use IdentityEvaluator<RT, 6>.
 *
 * \tparam RE Type of the underlying resultevalutor
 */
template <typename RE, typename MAT>
struct VanishingMaterialsWrapper
{
  using Underlying = RE;

  template <typename RE_>
  VanishingMaterialsWrapper(RE_&& resultEvaluator, const MAT& mat)
      : underlying_(std::forward<RE>(resultEvaluator)),
        mat_(mat) {}

  /**
   * \brief Calculate the result quantity by calculating the missing stress in z-direction and forwarding the result to
   * a resultevalutor
   * \param resultArray EigenVector containing the stress state in Voigt notation
   * \param comp component of result
   * \tparam R Type of the matrix
   * \return stress quantity
   */
  template <typename R>
  double operator()(const R& resultArray, const int comp) const {
    if constexpr (traits::isSpecializationNonTypeAndTypes<VanishingStrain, MAT>::value) {
      static_assert(R::CompileTimeTraits::RowsAtCompileTime == 3, "VanishingMaterialsWrapper is only valid for 2D.");
      auto nu                       = convertLameConstants(mat_.materialParameters()).toPoissonsRatio();
      auto sigZ                     = nu * (resultArray[0] + resultArray[1]);
      auto enlargedResultArray      = Eigen::Vector<double, 6>::Zero().eval();
      enlargedResultArray.head<2>() = resultArray.template head<2>();
      enlargedResultArray[2]        = sigZ;
      enlargedResultArray[5]        = resultArray[2];

      return underlying_(enlargedResultArray, comp);
    } else
      static_assert(Dune::AlwaysFalse<MAT>::value, "Material Type is not supported.");
  }
  /**
   * \brief Get the name of the result type
   * \return String representing the name
   */
  constexpr static std::string name() { return Underlying::name(); }

  /**
   * \brief Get the number of components in the result
   * \return Number of components
   */
  constexpr static int ncomps() { return Underlying::ncomps(); }

private:
  Underlying underlying_;
  MAT mat_;
};

template <typename ResultEvaluator, typename MAT>
VanishingMaterialsWrapper(ResultEvaluator&&, MAT&&) -> VanishingMaterialsWrapper<ResultEvaluator, MAT>;

/**
 * \brief Identity resultevalutor. Returns the results as is. Can for example be used with PlaneStrainWrapper.
 *
 * \tparam RT the requested result type
 * \tparam ncomps_ the amount of results in resultArray
 */
template <template <typename, int, int> class RT, int ncomps_>
struct IdentityEvaluator
{
  template <typename R>
  double operator()(const R& resultArray, const int comp) const {
    static_assert(R::CompileTimeTraits::RowsAtCompileTime >= ncomps_,
                  "Components in resultArray have to be at least as many as ncomps_.");
    return resultArray[comp];
  }

  constexpr static int ncomps() { return ncomps_; }

  constexpr static std::string name() { return toString<RT>(); }
};

/**
 * \brief Struct for calculating the 2d polar stress. This assumes cubed geometry.
 * \ingroup resultevaluators
 */
struct PolarStress
{
  /**
   * \brief Calculate the result quantity (von Mises stress)
   * \param resultArray EigenVector containing the stress state in Voigt notation
   * \param comp component of result
   * \tparam R Type of the matrix
   * \return von Mises stress
   */
  template <typename R>
  double operator()(const R& resultArray, const auto& pos, const int comp) const {
    static_assert(R::CompileTimeTraits::RowsAtCompileTime == 3, "PolarStress is only valid for 2D.");

    // Offset to center the coordinate system in the reference geometry
    auto center = Dune::ReferenceElements<double, 2>::cube().geometry<0>(0).center();
    auto theta  = std::atan2(pos[1] - center[1], pos[0] - center[0]);

    const auto s_x  = resultArray[0];
    const auto s_y  = resultArray[1];
    const auto s_xy = resultArray[2];

    auto hpl   = 0.5 * (s_x + s_y);
    auto hmi   = 0.5 * (s_x - s_y);
    auto cos2t = std::cos(2 * theta);
    auto sin2t = std::sin(2 * theta);

    Dune::FieldVector<double, 3> sigmaPolar;
    sigmaPolar[0] = hpl + hmi * cos2t + s_xy * sin2t;
    sigmaPolar[1] = hpl - hmi * cos2t - s_xy * sin2t;
    sigmaPolar[2] = -hmi * sin2t + s_xy * cos2t;

    return sigmaPolar[comp];
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
};

} // namespace Ikarus::ResultEvaluators
