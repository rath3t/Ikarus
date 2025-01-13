// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file vanishingstress.hh
 * \brief Defines the VanishingStress material model and related functions.
 * \ingroup  materials
 */

#pragma once

#include "vanishinghelpers.hh"

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverfactory.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/nonlinearoperator.hh>

namespace Ikarus {

/**
 * \brief VanishingStress material model that enforces stress components to be zero.
 * \ingroup materials
 * \tparam stressIndexPair An array of MatrixIndexPair representing fixed stress components.
 * \tparam MI The underlying material model.
 */
template <auto stressIndexPair, typename MI>
struct VanishingStress : public Material<VanishingStress<stressIndexPair, MI>>
{
  using Underlying         = MI; ///< The underlying material type.
  using MaterialParameters = typename Underlying::MaterialParameters;

  static constexpr auto fixedPairs        = stressIndexPair;                     ///< Array of fixed stress components.
  static constexpr auto freeVoigtIndices  = createfreeVoigtIndices(fixedPairs);  ///< Free Voigt indices.
  static constexpr auto fixedVoigtIndices = createFixedVoigtIndices(fixedPairs); ///< Fixed Voigt indices.
  static constexpr auto fixedDiagonalVoigtIndicesSize =
      countDiagonalIndices(fixedPairs);                                ///< Number of fixed diagonal indices.
  static constexpr auto freeStrains = freeVoigtIndices.size();         ///< Number of free strains.
  using ScalarType                  = typename Underlying::ScalarType; ///< Scalar type.
  static constexpr bool isAutoDiff  = Concepts::AutodiffScalar<ScalarType>;

  static constexpr auto strainTag              = Underlying::strainTag;            ///< Strain tag.
  static constexpr auto stressTag              = Underlying::stressTag;            ///< Stress tag.
  static constexpr auto tangentModuliTag       = Underlying::tangentModuliTag;     ///< Tangent moduli tag.
  static constexpr bool energyAcceptsVoigt     = Underlying::energyAcceptsVoigt;   ///< Energy accepts Voigt notation.
  static constexpr bool stressToVoigt          = true;                             ///< Stress to Voigt notation.
  static constexpr bool stressAcceptsVoigt     = true;                             ///< Stress accepts Voigt notation.
  static constexpr bool moduliToVoigt          = true;                             ///< Moduli to Voigt notation.
  static constexpr bool moduliAcceptsVoigt     = true;                             ///< Moduli accepts Voigt notation.
  static constexpr double derivativeFactorImpl = Underlying::derivativeFactorImpl; ///< Derivative factor.

  /**
   * \brief Constructor for VanishingStress.
   * \param mat The underlying material model.
   * \param tol Tolerance for stress reduction.
   */
  explicit VanishingStress(MI mat, typename MI::ScalarType tol = 1e-12)
      : matImpl_{mat},
        tol_{tol} {}

  [[nodiscard]] constexpr static std::string nameImpl() noexcept {
    auto matName = MI::name() + "_VanishingStress(";
    for (auto p : fixedPairs)
      matName += "(" + std::to_string(p.row) + std::to_string(p.col) + ")";
    matName += ")";
    return matName;
  }

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return matImpl_.materialParametersImpl(); }

  /**
   * \brief Computes the stored energy for the VanishingStress material.
   * \tparam Derived The derived type of the input matrix.
   * \param E The Green-Lagrangian strain.
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& E) const {
    const auto [nonOp, Esol] = reduceStress(E);
    return matImpl_.storedEnergyImpl(Esol);
  }

  /**
   * \brief Computes the stresses for the VanishingStress material.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The Green-Lagrangian strain.
   * \return StressMatrix The stresses.
   */
  template <bool voigt, typename Derived>
  auto stressesImpl(const Eigen::MatrixBase<Derived>& E) const {
    const auto [nonOp, Esol] = reduceStress(E);
    auto stressesRed         = matImpl_.template stresses<Underlying::strainTag, true>(Esol);

    if constexpr (voigt) {
      return removeCol(stressesRed, fixedVoigtIndices);
    } else
      return fromVoigt(stressesRed, false);
  }

  /**
   * \brief Computes the tangent moduli for the VanishingStress material.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The Green-Lagrangian strain.
   * \return TangentModuli The tangent moduli.
   */
  template <bool voigt, typename Derived>
  auto tangentModuliImpl(const Eigen::MatrixBase<Derived>& E) const {
    const auto [nonOp, Esol] = reduceStress(E);
    auto C                   = matImpl_.template tangentModuli<Underlying::strainTag, true>(Esol);
    if constexpr (voigt)
      return staticCondensation(C, fixedVoigtIndices);
    else
      return fromVoigt(C);
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam ScalarTypeOther The target scalar type.
   * \return VanishingStress The rebound VanishingStress material.
   */
  template <typename ScalarTypeOther>
  auto rebind() const {
    auto reboundMatImpl = matImpl_.template rebind<ScalarTypeOther>();
    return VanishingStress<stressIndexPair, decltype(reboundMatImpl)>(reboundMatImpl, tol_);
  }

  /**
   * \brief Returns a const reference to the underlying material.
   */
  auto& underlying() const { return matImpl_; }

private:
  /**
   * \brief Initializes unknown strains based on fixed indices.
   * \tparam Derived The derived type of the input matrix.
   * \param E The input strain matrix.
   */
  template <typename Derived>
  void initUnknownStrains(Eigen::MatrixBase<Derived>& E) const {
    for (size_t i = 0; i < fixedPairs.size(); ++i) {
      ScalarType initialVal = E(fixedPairs[i].row, fixedPairs[i].col);
      if constexpr (strainTag == StrainTags::deformationGradient or strainTag == StrainTags::rightCauchyGreenTensor) {
        if (Dune::FloatCmp::eq(initialVal, ScalarType(0.0)) and (fixedPairs[i].row == fixedPairs[i].col))
          initialVal = ScalarType(1.0);
      }
      if (fixedPairs[i].row != fixedPairs[i].col)
        initialVal = ScalarType(0.0);
      E(fixedPairs[i].row, fixedPairs[i].col) = E(fixedPairs[i].col, fixedPairs[i].row) = initialVal;
    }
  }

  /**
   * \brief Reduces stress components to satisfy the vanishing stress condition.
   * \tparam Derived The derived type of the input matrix.
   * \param Eraw The input strain matrix.
   * \return std::pair<NonLinearOperator, decltype(auto)> The stress reduction result.
   */
  template <typename Derived>
  auto reduceStress(const Eigen::MatrixBase<Derived>& Eraw) const {
    auto E = Impl::maybeFromVoigt(Eraw);
    initUnknownStrains(E);

    std::array<size_t, fixedDiagonalVoigtIndicesSize> fixedDiagonalVoigtIndices;
    for (size_t ri = 0; auto i : fixedVoigtIndices) {
      auto indexPair = fromVoigt(i);
      if (indexPair[0] == indexPair[1])
        fixedDiagonalVoigtIndices[ri++] = i;
    }

    auto f = [&](auto&) {
      auto S = matImpl_.template stresses<Underlying::strainTag, true>(E);
      return S(fixedDiagonalVoigtIndices).eval();
    };
    auto df = [&](auto&) {
      auto moduli = (matImpl_.template tangentModuli<Underlying::strainTag, true>(E)).eval();
      return (moduli(fixedDiagonalVoigtIndices, fixedDiagonalVoigtIndices) / Underlying::derivativeFactor).eval();
    };

    auto Er    = E(fixedDiagonalVoigtIndices, fixedDiagonalVoigtIndices).eval().template cast<ScalarType>();
    auto nonOp = Ikarus::NonLinearOperator(functions(f, df), parameter(Er));

    auto linearSolver   = [](auto& r, auto& A) { return (A.inverse() * r).eval(); };
    auto updateFunction = [&](auto& /* Ex33 */, auto& ecomps) {
      for (int ri = 0; auto i : fixedDiagonalVoigtIndices) {
        auto indexPair = fromVoigt(i);
        E(indexPair[0], indexPair[1]) += ecomps(ri++);
      }
    };

    int minIter = isAutoDiff ? 1 : 0;
    // THE CTAD is broken for designated initializers in clang 16, when we drop support this can be simplified
    NewtonRaphsonConfig<decltype(linearSolver), decltype(updateFunction)> nrs{
        .parameters     = {.tol = tol_, .maxIter = 100, .minIter = minIter},
        .linearSolver   = linearSolver,
        .updateFunction = updateFunction
    };

    auto nr = createNonlinearSolver(std::move(nrs), nonOp);
    if (!static_cast<bool>(nr->solve()))
      DUNE_THROW(Dune::MathError, "The stress reduction of material " << nameImpl() << " was unsuccessful\n"
                                                                      << "The strains are\n"
                                                                      << E << "\n The stresses are\n"
                                                                      << f(Er));
    return std::make_pair(nonOp, E);
  }

  Underlying matImpl_; ///< The underlying material model.
  double tol_{};       ///< Tolerance for stress reduction.
};

/**
 * \brief Factory function to create a VanishingStress material with specified stress indices.
 * \tparam matrixIndexPair The array of MatrixIndexPair representing fixed stress components.
 * \tparam MaterialImpl The underlying material model.
 * \param mat The underlying material model.
 * \param p_tol Tolerance for stress reduction.
 * \return VanishingStress The created VanishingStress material.
 */
template <Impl::MatrixIndexPair... stressIndexPair, typename MaterialImpl>
auto makeVanishingStress(MaterialImpl mat, typename MaterialImpl::ScalarType p_tol = 1e-12) {
  return VanishingStress<std::to_array({stressIndexPair...}), MaterialImpl>(mat, p_tol);
}

/**
 * \brief Factory function to create a VanishingStress material for plane stress conditions.
 * \tparam MaterialImpl The underlying material model.
 * \param mat The underlying material model.
 * \param tol Tolerance for stress reduction.
 * \return VanishingStress The created VanishingStress material for plane stress.
 */
template <typename MaterialImpl>
auto planeStress(const MaterialImpl& mat, typename MaterialImpl::ScalarType tol = 1e-8) {
  return makeVanishingStress<Impl::MatrixIndexPair{2, 1}, Impl::MatrixIndexPair{2, 0}, Impl::MatrixIndexPair{2, 2}>(
      mat, tol);
}

/**
 * \brief Factory function to create a VanishingStress material for a shell material with zero normal stress
 * condition.
 * \tparam MaterialImpl The underlying material model.
 * \param mat The underlying material model.
 * \param tol Tolerance for stress reduction.
 * \return VanishingStress The created VanishingStress material for plane stress.
 */
template <typename MaterialImpl>
auto shellMaterial(const MaterialImpl& mat, typename MaterialImpl::ScalarType tol = 1e-8) {
  return makeVanishingStress<Impl::MatrixIndexPair{2, 2}>(mat, tol);
}

/**
 * \brief Factory function to create a VanishingStress material for a beam material with two zero normal stress
 * condition.
 * \tparam MaterialImpl The underlying material model.
 * \param mat The underlying material model.
 * \param p_tol Tolerance for stress reduction.
 * \return VanishingStress The created VanishingStress material for plane stress.
 */
template <typename MaterialImpl>
auto beamMaterial(const MaterialImpl& mat, typename MaterialImpl::ScalarType tol = 1e-8) {
  return makeVanishingStress<Impl::MatrixIndexPair{1, 1}, Impl::MatrixIndexPair{2, 2}>(mat, tol);
}
} // namespace Ikarus
