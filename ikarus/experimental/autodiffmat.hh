// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file autodiffmat.hh
 * \brief Implementation of the AutoDiff-based material model.
 * \ingroup materials
 */

#pragma once

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <autodiff/forward/utils/derivative.hpp>

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Experimental {

/**
 * \brief Implementation of a AutoDiff-based material model.
 * \ingroup experimental
 * \details As of now using autodiff on the implemented materials doesn't always yield the expected result. Therefore
 * this is labelled as an experimental feature for now. Please validate your results independently.
 *
 * \tparam RealMAT Type of the original material model.
 * \tparam forceAutoDiffV Force automatic differentiation to compute tangentModuli via stresses and stresses via energy.
 * \tparam forceAutoDiffS Force automatic differentiation to compute tangentModuli and stresses via energy.
 */
template <typename RealMAT, bool forceAutoDiffV = false, bool forceAutoDiffS = false>
struct AutoDiffMAT : public RealMAT
{
  using ScalarType     = typename RealMAT::ScalarType;
  using StrainMatrix   = typename RealMAT::StrainMatrix;
  using StressMatrix   = typename RealMAT::StressMatrix;
  using MaterialTensor = typename RealMAT::MaterialTensor;

  using MaterialParameters = typename RealMAT::MaterialParameters;

  static constexpr int dim           = RealMAT::dim;
  static constexpr int nVoigtIndices = dim * (dim + 1) / 2;

  static constexpr auto strainTag              = RealMAT::strainTag;
  static constexpr auto stressTag              = RealMAT::stressTag;
  static constexpr auto tangentModuliTag       = RealMAT::tangentModuliTag;
  static constexpr bool energyAcceptsVoigt     = RealMAT::energyAcceptsVoigt;
  static constexpr bool stressToVoigt          = RealMAT::stressToVoigt;
  static constexpr bool stressAcceptsVoigt     = RealMAT::stressAcceptsVoigt;
  static constexpr bool moduliToVoigt          = RealMAT::moduliToVoigt;
  static constexpr bool moduliAcceptsVoigt     = RealMAT::moduliAcceptsVoigt;
  static constexpr double derivativeFactorImpl = RealMAT::derivativeFactorImpl;

  /**
   * \brief Constructor for the AutoDiffMAT class.
   * Forward the construction to the underlying element
   *
   * \tparam Args Variadic template for constructor arguments.
   * \param args Constructor arguments.
   */
  template <typename... Args>
  explicit AutoDiffMAT(Args&&... args)
      : RealMAT{std::forward<Args>(args)...} {}

  [[nodiscard]] constexpr static std::string nameImpl() noexcept { return "AutoDiff: " + RealMAT::name(); }

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return realMAT().materialParametersImpl(); }

  /**
   * \brief Computes the stored energy in the underlying material model.
   * \tparam Derived The derived type of the input matrix.
   * \param E The strain mesasure.
   * \return The stored energy.
   */
  template <StrainTags tag, typename Derived>
  auto storedEnergy(const Eigen::MatrixBase<Derived>& E) const {
    if constexpr (requires { realMAT().template storedEnergy<tag>(E); }) {
      auto mat_ad = realMAT().template rebind<autodiff::dual>();
      return mat_ad.template storedEnergy<tag>(E);
    } else {
      static_assert(Dune::AlwaysFalse<AutoDiffMAT>::value,
                    "Appropriate storedEnergy function not is implemented for the chosen material model.");
    }
  }

  /**
   * \brief Computes the stresses in the underlying material model.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The strain mesasure.
   * \return The stresses.
   */
  template <StrainTags tag, bool voigt = true, typename Derived>
  auto stresses(const Eigen::MatrixBase<Derived>& E) const {
    if constexpr (requires { realMAT().template stresses<tag>(E); } and not(forceAutoDiffV or forceAutoDiffS)) {
      return realMAT().template stresses<tag>(E);
    } else if constexpr (requires { realMAT().template storedEnergy<tag>(E); }) {
      static_assert(!Concepts::EigenVector<Derived>,
                    "The strain measure used for autodiff has to be in matrix notation.");
      auto mat_ad = realMAT().template rebind<autodiff::dual>();

      auto f = [&](const auto& x) { return mat_ad.template storedEnergy<tag>(x); };

      auto dx = Eigen::Vector<autodiff::dual, nVoigtIndices>{};

      dx = toVoigt(E.derived());
      autodiff::dual e;

      auto g = Eigen::Vector<double, nVoigtIndices>{};
      gradient(f, autodiff::wrt(dx), autodiff::at(dx), e, g);

      return (derivativeFactorImpl * g).eval();
    } else {
      static_assert(Dune::AlwaysFalse<AutoDiffMAT>::value,
                    "Appropriate storedEnergy function not is implemented for the chosen material model.");
    }
  }

  /**
   * \brief Computes the tangent moduli in the underlying material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The strain mesasure.
   * \return The tangent moduli.
   */
  template <StrainTags tag, bool voigt = true, typename Derived>
  auto tangentModuli(const Eigen::MatrixBase<Derived>& E) const {
    if constexpr (requires { realMAT().template tangentModuli<tag>(E); } and not(forceAutoDiffV or forceAutoDiffS)) {
      return realMAT().template tangentModuli<tag>(E);
    } else if constexpr (requires { realMAT().template stresses<tag>(E); } and forceAutoDiffV and not forceAutoDiffS) {
      static_assert(!Concepts::EigenVector<Derived>,
                    "The strain measure used for autodiff has to be in matrix notation.");
      auto mat_ad = realMAT().template rebind<autodiff::dual>();

      auto f = [&](const auto& x) { return mat_ad.template stresses<tag>(x); };

      auto dx = Eigen::Vector<autodiff::dual, nVoigtIndices>{};

      dx = toVoigt(E.derived());
      Eigen::VectorXdual g(nVoigtIndices);

      auto h = Eigen::Matrix<double, nVoigtIndices, nVoigtIndices>{};
      jacobian(f, autodiff::wrt(dx), autodiff::at(dx), g, h);

      return (derivativeFactorImpl * h).eval();
    } else if constexpr (requires { realMAT().template storedEnergy<tag>(E); }) {
      static_assert(!Concepts::EigenVector<Derived>,
                    "The strain measure used for autodiff has to be in matrix notation.");
      auto mat_ad = realMAT().template rebind<autodiff::dual2nd>();

      auto f = [&](const auto& x) { return mat_ad.template storedEnergy<tag>(x); };

      Eigen::Matrix<autodiff::dual2nd, nVoigtIndices, 1> dx = toVoigt(E.derived());

      autodiff::dual2nd e;
      Eigen::Matrix<double, nVoigtIndices, 1> g;
      Eigen::Matrix<double, nVoigtIndices, nVoigtIndices> h;

      h = autodiff::hessian(f, autodiff::wrt(dx), autodiff::at(dx), e, g);

      return (derivativeFactorImpl * derivativeFactorImpl * h).eval();
    } else {
      static_assert(
          Dune::AlwaysFalse<AutoDiffMAT>::value,
          "Appropriate storedEnergy and stresses functions are not implemented for the chosen material model.");
    }
  }

  /**
   * \brief Get the reference to the base material.
   *
   * \return The reference to the base material.
   */
  const RealMAT& realMAT() const { return *this; }

private:
  MaterialParameters materialParameter_;
};

} // namespace Ikarus::Experimental
