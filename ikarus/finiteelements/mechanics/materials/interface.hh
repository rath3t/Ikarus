// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file interface.hh
 * \brief Contains the Material interface class and related template functions for material properties.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/strainconversions.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/linearalgebrahelper.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus::Materials {

#ifndef DOXYGEN
template <class MImpl>
struct Material;

template <auto stressIndexPair, typename MImpl>
struct VanishingStress;

template <auto strainIndexPair, typename MImpl>
struct VanishingStrain;
#endif

/**
 * \brief Template function for checking if the strain size is correct.
 *
 * The given strain quantity has to be a Eigen::Vector6 or a Eigen::Matrix3
 *
 * \tparam MAT Type of the material.
 * \tparam S Type of the strains.
 */
template <typename MAT, typename S>
consteval bool hasCorrectSize() {
  if constexpr (Concepts::EigenVector6<S> or Concepts::EigenMatrix33<S>)
    return true;
  if constexpr (MAT::isReduced and Concepts::EigenVector<S>) {
    return S::RowsAtCompileTime == MAT::freeStrains;
  } else
    return false;
}

/**
 * \brief Template concept for ensuring correct strain size.
 *
 * \tparam MAT Type of the material.
 * \tparam S Type of the strains.
 */
template <typename MAT, typename S>
concept CorrectStrainSize = hasCorrectSize<MAT, S>();

/**
 * \brief Interface classf or materials.
 * \ingroup materials
 *    \details Consider a hyper elastic material with the free Helmholtz energy
 *        \f[\require{cases}\psi: \begin{cases}\mathbb{R}^{3\times 3} &\rightarrow \mathbb{R} \\ \BC
 * &\mapsto \psi(\BC) \end{cases}.\f]
 *
 * Then, the value of this potential energy is return by \link Material< MaterialImpl >::storedEnergy storedEnergy
 * \endlink and its first derivative (the stresses) by \link Material< MaterialImpl >::stresses stresses \endlink and
 * the second derivatives (the tangent moduli) by \link Material< MaterialImpl >::tangentModuli tangentModuli
 * \endlink.
 *
 * The passed strains can be in several formats, i.e.,
 *   \f$\BC\f$ can be the [right Cauchy-Green
 * tensor](https://en.wikipedia.org/wiki/Finite_strain_theory#Cauchy_strain_tensor_(right_Cauchy%E2%80%93Green_deformation_tensor)),
 * the [deformation gradient](https://en.wikipedia.org/wiki/Finite_strain_theory#Deformation_gradient_tensor)
 * \f$\mathbf{F}\f$ or linear strains. The current supported tags are given by Ikarus::StrainTags.
 * \tparam MI Type of the underlying material implementation.
 */
template <class MI>
struct Material
{
  using MaterialImpl = MI; ///< Type of material implementation

  /**
   * \brief Static constant for determining if the material has vanishing stress or strain components (is reduced).
   */
  static constexpr bool isReduced = traits::isSpecializationNonTypeAndTypes<VanishingStress, MaterialImpl>::value or
                                    traits::isSpecializationNonTypeAndTypes<VanishingStrain, MaterialImpl>::value;

  /**
   * \brief Const accessor to the underlying material (CRTP).
   *
   * \return Const reference to the underlying material.
   */
  constexpr const MaterialImpl& impl() const { return static_cast<const MaterialImpl&>(*this); }

  /**
   * \brief Accessor to the underlying material (CRTP).
   *
   * \return Reference to the underlying material.
   */
  constexpr MaterialImpl& impl() { return static_cast<MaterialImpl&>(*this); }

  /**
   * \brief Get the name of the implemented material.
   *
   * \return Name of the material.
   */
  [[nodiscard]] constexpr static std::string name() { return MI::nameImpl(); }

  /**
   * \brief Returns the material parameters stored in the implemented material.
   * \return Material parameter.
   */
  [[nodiscard]] auto materialParameters() const { return impl().materialParametersImpl(); }

  /**
   * \brief This factor denotes the differences between the returned stresses and moduli and the passed strain
   * \details For neoHooke the inserted quantity is $C$ the Green-Lagrangian strain tensor, the function relation
   * between the energy and the stresses is $S = 1 \dfrac{\partial\psi(E)}{\partial E}$. This factor is the pre factor,
   * which is the difference to the actual derivative, which is written here
   */
  static constexpr double derivativeFactor = MI::derivativeFactorImpl;

  /**
   * \brief Return the stored potential energy of the material.
   *
   *\details This function return the free Helmholtz energy of the material
   *
   * \tparam tag Strain tag indicating which strain tensor components are passed.
   * \tparam Derived The underlying Eigen type.
   * \param Eraw The strain tensor components passed in Voigt notation or matrix notation.
   * \return Scalar return of stored energy.
   */
  template <StrainTags tag, typename Derived>
  requires CorrectStrainSize<MaterialImpl, Derived>
  [[nodiscard]] auto storedEnergy(const Eigen::MatrixBase<Derived>& Eraw) const {
    decltype(auto) Ev = enlargeIfReduced<Material>(Eraw);
    decltype(auto) E  = transformStrain<tag, MaterialImpl::strainTag>(Ev);

    if constexpr (Concepts::EigenVector<Derived>) { // receiving vector means voigt notation
      if constexpr (MaterialImpl::energyAcceptsVoigt)
        return impl().storedEnergyImpl(toVoigt(E));
      else
        return impl().storedEnergyImpl(E);
    } else
      return impl().storedEnergyImpl(E);
  }

  /**
   * \brief Get the stresses of the material.
   *
   * \tparam tag Strain tag indicating which strain tensor components are passed.
   * \tparam voigt Boolean indicating whether to return Voigt-shaped result.
   * \tparam Derived The underlying Eigen type.
   * \param Eraw The strain tensor components passed in Voigt notation or matrix notation.
   * \return Vectorial or Matrix return of stresses.
   */
  template <StrainTags tag, bool voigt = true, typename Derived>
  requires CorrectStrainSize<MaterialImpl, Derived>
  [[nodiscard]] auto stresses(const Eigen::MatrixBase<Derived>& Eraw) const {
    decltype(auto) Ev = enlargeIfReduced<Material>(Eraw);
    decltype(auto) E  = transformStrain<tag, MaterialImpl::strainTag>(Ev);
    if constexpr (voigt and MaterialImpl::stressToVoigt == false)
      // user requests a Voigt shaped return but material is not able to. Therefore, we transform it here.
      return toVoigt(stressesMaybeTransformInputToVoigt<false>(E), false);
    else
      return stressesMaybeTransformInputToVoigt<voigt>(E);
  }

  /**
   * \brief Get the tangentModuli of the material.
   *
   * \tparam tag Strain tag indicating which strain tensor components are passed.
   * \tparam voigt Boolean indicating whether to return Voigt-shaped result.
   * \tparam Derived The underlying Eigen type.
   * \param Eraw The strain tensor components passed in Voigt notation or matrix notation.
   * \return Tangent moduli in Voigt notation or as fourth-order tensor.
   */
  template <StrainTags tag, bool voigt = true, typename Derived>
  requires CorrectStrainSize<MaterialImpl, Derived>
  [[nodiscard]] auto tangentModuli(const Eigen::MatrixBase<Derived>& Eraw) const {
    decltype(auto) Ev = enlargeIfReduced<Material>(Eraw);
    decltype(auto) E  = transformStrain<tag, MaterialImpl::strainTag>(Ev);
    if constexpr (voigt and MaterialImpl::moduliToVoigt == false)
      // user request a Voigt shaped return but material is not able to. Therefore, we transform it here.
      return toVoigt(tangentModuliMaybeTransformInputToVoigt<false>(E));
    else
      return tangentModuliMaybeTransformInputToVoigt<voigt>(E);
  }

  /**
   * \brief Rebind material to a different scalar type.
   *
   * Useful for using automatic differentiation.
   *
   * \tparam STO The scalar type to rebind to.
   * \return Rebound material.
   */
  template <typename STO>
  auto rebind() const {
    return impl().template rebind<STO>();
  }

private:
  template <bool voigt = true, typename Derived>
  auto stressesMaybeTransformInputToVoigt(const Eigen::MatrixBase<Derived>& E) const {
    if constexpr (Concepts::EigenVector<Derived>) { // receiving vector means Voigt notation
      if constexpr (MaterialImpl::stressAcceptsVoigt)
        return impl().template stressesImpl<voigt>(E);
      else // material is not able to accept Voigt shaped Input. Therefore, we transform it before.
        return impl().template stressesImpl<voigt>(fromVoigt(E.derived()));
    } else
      return impl().template stressesImpl<voigt>(E.derived());
  }

  template <bool voigt = true, typename Derived>
  auto tangentModuliMaybeTransformInputToVoigt(const Eigen::MatrixBase<Derived>& E) const {
    if constexpr (Concepts::EigenVector<Derived>) { // receiving vector means voigt notation
      if constexpr (MaterialImpl::moduliAcceptsVoigt)
        return impl().template tangentModuliImpl<voigt>(E);
      else
        return impl().template tangentModuliImpl<voigt>(fromVoigt(E.derived()));
    } else
      return impl().template tangentModuliImpl<voigt>(E.derived());
  }
};

} // namespace Ikarus::Materials
