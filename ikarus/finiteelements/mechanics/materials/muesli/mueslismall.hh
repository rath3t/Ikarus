// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file Muesli.hh
 * \brief Implementation of the Muesli material model.
 * \ingroup materials
 */

#pragma once

#include <muesli/muesli.h>

#include <Eigen/Eigen>

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/muesli/mueslihelpers.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials::Muesli {

/**
 * \brief Wrapper class for small strain materials from the muesli library. It can be templated with all materials
 * derived from muesli::elasticIsotropicMaterial. It models the Ikarus material interface.
 *
 * \tparam SM muesli material model implementation
 */
template <typename SM = muesli::elasticIsotropicMaterial>
requires(std::is_base_of_v<muesli::smallStrainMaterial, SM>)
struct SmallStrain : public Material<SmallStrain<SM>>
{
  using MaterialModel                 = SM;
  using ScalarType                    = double;
  static constexpr int worldDimension = 3;
  using StrainMatrix                  = Eigen::Matrix<ScalarType, worldDimension, worldDimension>;
  using StressMatrix                  = StrainMatrix;
  using MaterialParameters            = Muesli::MaterialProperties;

  static constexpr auto strainTag              = StrainTags::linear;
  static constexpr auto stressTag              = StressTags::linear;
  static constexpr auto tangentModuliTag       = TangentModuliTags::Material;
  static constexpr bool energyAcceptsVoigt     = false;
  static constexpr bool stressToVoigt          = false;
  static constexpr bool stressAcceptsVoigt     = false;
  static constexpr bool moduliToVoigt          = false;
  static constexpr bool moduliAcceptsVoigt     = false;
  static constexpr double derivativeFactorImpl = 1;

  [[nodiscard]] constexpr static std::string nameImpl() noexcept { return "Muesli_SmallStrain: " + materialName<SM>(); }

  /**
   * \brief Constructor for SmallStrain muesli materials (only activated for isotropic linear elasticity).
   * \param mpt Arbitrary Material Parameter tuple defined in physicshelper.hh
   */
  template <Concepts::MPTuple MPT>
  requires(std::same_as<MaterialModel, muesli::elasticIsotropicMaterial>)
  explicit SmallStrain(const MPT& mpt)
      : materialParameter_{Muesli::propertiesFromIkarusMaterialParameters(mpt)},
        material_{Dune::className<SM>(), materialParameter_},
        mp_{material_.createMaterialPoint()} {}

  /**
   * \brief CConstructor for SmallStrain muesli materials
   * \param mpt Muesli materialproperties
   */
  explicit SmallStrain(const MaterialParameters& mpt)
      : materialParameter_{mpt},
        material_{Dune::className<SM>(), mpt},
        mp_{material_.createMaterialPoint()} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameter_; }

  /**
   * \brief Computes the stored energy in the Neo-Hookean material model.
   * \tparam Derived The derived type of the input matrix.
   * \param E The linear strain tensor.
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& E) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!Concepts::EigenVector<Derived>) {
      updateState(E);
      return mp_->storedEnergy();
    } else
      static_assert(!Concepts::EigenVector<Derived>,
                    "MuesliSmallStrain energy can only be called with a matrix and not a vector in Voigt notation");
  }

  /**
   * \brief Computes the stresses in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The linear strain tensor.
   * \return StressMatrix The stresses.
   */
  template <bool voigt, typename Derived>
  auto stressesImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      if constexpr (!Concepts::EigenVector<Derived>) {
        updateState(C);
        mp_->stress(stress_);
        return Muesli::toMatrix(stress_);
      } else
        static_assert(!Concepts::EigenVector<Derived>,
                      "MuesliSmallStrain can only be called with a matrix and not a vector in Voigt notation");
    } else
      static_assert(voigt == false, "MuesliSmallStrain does not support returning stresses in Voigt notation");
  }

  /**
   * \brief Computes the tangent moduli in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The linear strain tensor.
   * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> The tangent moduli.
   */
  template <bool voigt, typename Derived>
  auto tangentModuliImpl(const Eigen::MatrixBase<Derived>& E) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      if constexpr (!Concepts::EigenVector<Derived>) {
        updateState(E);
        mp_->tangentTensor(tangentModuli_);
        return Muesli::toTensor(tangentModuli_);
      } else
        static_assert(!Concepts::EigenVector<Derived>,
                      "MuesliSmallStrain can only be called with a matrix and not a vector in Voigt notation");
    } else
      static_assert(voigt == false, "MuesliSmallStrain does not support returning tangent moduli in Voigt notation");
  }

  /**
   * \brief Returns the underlying muesli material implementation
   * \return auto& reference to the musli material
   */
  auto& material() const { return material_; }

  /**
   * \brief asserts that the materialpoint pointer is not null
   */
  bool assertMP() const { return mp_.get() != NULL; }

  SmallStrain(const SmallStrain& other)
      : materialParameter_{other.materialParameter_},
        material_{Dune::className<SM>(), materialParameter_},
        mp_{material_.createMaterialPoint()} {}

private:
  MaterialParameters materialParameter_;
  MaterialModel material_;
  std::unique_ptr<muesli::smallStrainMP> mp_;

  mutable istensor strain_{};
  mutable istensor stress_{};
  mutable itensor4 tangentModuli_{};

  template <typename Derived>
  void updateState(const Eigen::MatrixBase<Derived>& E) const {
    Muesli::toistensor(strain_, E);
    mp_->updateCurrentState(0.0, strain_);
  }
};

} // namespace Ikarus::Materials::Muesli
