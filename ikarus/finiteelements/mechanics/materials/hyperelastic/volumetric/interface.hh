// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file interface.hh
 * \brief Implementation of the volumetric part of a hyperelastic material.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/concepts.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Interface for the volumetric part of a hyperelastic material. Has to be parametrized with a volumetric
 * function.
 * \ingroup materials
 *
 * \details The volumetric part of the hyperelastic model, i.e., related to \f$ U(J) \f$, is
 * parametrized with a certain volumetric function (VF) implemented in terms of the determinant of the deformation
 * gradient. The underlying volumetric function must only implement the energy \f$ U(J) \f$ and its first
 * and second derivatives w.r.t \f$ J \f$.
 *
 * \tparam VF volumetric function, has to adhere to the \ref VolumetricFunction.
 */
template <Concepts::VolumetricFunction VF>
struct Volumetric
{
  using ScalarType = typename VF::ScalarType;
  using JType      = typename VF::JType;

  using VolumetricFunction = VF;
  using MaterialParameter  = double;

  [[nodiscard]] constexpr static std::string name() noexcept { return "Volumetric function: " + VF::name(); }

  /**
   * \brief Construct a new volumetric part
   *
   * \param matPar Materialparameter that is considered as the penalty parameter for the constraint against volumetric
   * deformations. Typically, if the deviatoric part of the energy function is a function of the toal principal
   * stretches, then this material parameter is the  Lamé's first parameter and if the energy is a function of the
   * deviatoric principal stretches, then bulk modulus is used.
   *
   * \param vf the volumetric function
   */
  Volumetric(MaterialParameter matPar, const VF vf)
      : matPar_{matPar},
        volumetricFunction_{vf} {}

  /**
   * \brief Computes stored energy of the volumetric function
   *
   * \param J determinant of the deformation gradient \f$ J = \det\BF \f$
   * \return ScalarType energy
   */
  ScalarType storedEnergy(const JType& J) const { return matPar_ * volumetricFunction_.storedEnergyImpl(J); };

  /**
   * \brief Computes the first derivatives of the energy of the volumetric function w.r.t \f$ J \f$
   *
   * \param J determinant of the deformation gradient \f$ J = \det\BF \f$
   * \return ScalarType energy
   */
  ScalarType firstDerivative(const JType& J) const { return matPar_ * volumetricFunction_.firstDerivativeImpl(J); }

  /**
   * \brief Computes the second derivatives of the energy of the volumetric function w.r.t \f$ J \f$
   *
   * \param J determinant of the deformation gradient \f$ J = \det\BF \f$
   * \return ScalarType energy
   */
  ScalarType secondDerivative(const JType& J) const { return matPar_ * volumetricFunction_.secondDerivativeImpl(J); };

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return The rebound volumetric part.
   */
  template <typename STO>
  auto rebind() const {
    auto reboundVF = volumetricFunction_.template rebind<STO>();
    return Volumetric<decltype(reboundVF)>(matPar_, reboundVF);
  }

private:
  MaterialParameter matPar_;
  VF volumetricFunction_;
};
} // namespace Ikarus::Materials
