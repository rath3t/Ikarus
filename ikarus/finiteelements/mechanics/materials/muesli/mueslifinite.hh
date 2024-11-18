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

template <typename FM>
requires(std::is_base_of_v<muesli::finiteStrainMaterial, FM>)
struct FiniteStrain : public Material<FiniteStrain<FM>>
{
  using MaterialModel                 = FM;
  using ScalarType                    = double;
  static constexpr int worldDimension = 3;
  using StrainMatrix                  = Eigen::Matrix<ScalarType, worldDimension, worldDimension>;
  using StressMatrix                  = StrainMatrix;
  using MaterialParameters            = Muesli::MaterialProperties;

  static constexpr auto strainTag              = StrainTags::rightCauchyGreenTensor;
  static constexpr auto stressTag              = StressTags::PK2;
  static constexpr auto tangentModuliTag       = TangentModuliTags::Material;
  static constexpr bool energyAcceptsVoigt     = false;
  static constexpr bool stressToVoigt          = false;
  static constexpr bool stressAcceptsVoigt     = false;
  static constexpr bool moduliToVoigt          = false;
  static constexpr bool moduliAcceptsVoigt     = false;
  static constexpr double derivativeFactorImpl = 1;

  [[nodiscard]] constexpr static std::string nameImpl() noexcept { return "FiniteStrain: " + materialName<FM>(); }

  /**
   * \briefCConstructor for FiniteStrain muesli materials
   * \param mpt Muesli materialproperties
   */
  explicit FiniteStrain(const MaterialParameters& mpt)
      : materialParameter_{mpt},
        material_{Dune::className<FM>(), mpt},
        mp_{material_.createMaterialPoint()} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameter_; }

  /**
   * \brief Computes the stored energy in the Neo-Hookean material model.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!Concepts::EigenVector<Derived>) {
      updateState(C);
      return mp_->storedEnergy();

    } else
      static_assert(!Concepts::EigenVector<Derived>,
                    "MuesliFiniteStrain energy can only be called with a matrix and not a vector in Voigt notation");
  }

  /**
   * \brief Computes the stresses in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return StressMatrix The stresses.
   */
  template <bool voigt, typename Derived>
  auto stressesImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      if constexpr (!Concepts::EigenVector<Derived>) {
        updateState(C);
        mp_->secondPiolaKirchhoffStress(stress_);
        return Muesli::toMatrix<ScalarType>(stress_);
      } else
        static_assert(!Concepts::EigenVector<Derived>,
                      "MuesliFiniteStrain can only be called with a matrix and not a vector in Voigt notation");
    } else
      static_assert(voigt == false, "MuesliFiniteStrain does not support returning stresses in Voigt notation");
  }

  /**
   * \brief Computes the tangent moduli in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> The tangent moduli.
   */
  template <bool voigt, typename Derived>
  auto tangentModuliImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      if constexpr (!Concepts::EigenVector<Derived>) {
        updateState(C);

        mp_->convectedTangent(tangentModuli_);
        return Muesli::toTensor<ScalarType>(tangentModuli_);
      } else
        static_assert(!Concepts::EigenVector<Derived>,
                      "MuesliFiniteStrain can only be called with a matrix and not a vector in Voigt notation");
    } else
      static_assert(voigt == false, "MuesliFiniteStrain does not support returning tangent moduli in Voigt notation");
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

  FiniteStrain(const FiniteStrain& other)
      : materialParameter_{other.materialParameter_},
        material_{Dune::className<FM>(), materialParameter_},
        mp_{material_.createMaterialPoint()} {}

private:
  MaterialParameters materialParameter_;
  MaterialModel material_;
  std::unique_ptr<muesli::finiteStrainMP> mp_;

  mutable istensor strain_{};
  mutable istensor stress_{};
  mutable itensor4 tangentModuli_{};

  template <typename Derived>
  void updateState(const Eigen::MatrixBase<Derived>& C) const {
    Muesli::toistensor(strain_, transformStrain<strainTag, StrainTags::deformationGradient>(C));
    mp_->updateCurrentState(0.0, strain_);
  }
};

} // namespace Ikarus::Materials::Muesli
