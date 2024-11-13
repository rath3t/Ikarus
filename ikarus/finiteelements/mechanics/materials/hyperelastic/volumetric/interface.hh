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
 * \brief Interface for the volumetric part opf a hyperelastic material. Has to be parametrized with a volumetric
 * function.
 * \tparam VF volumetric function, has to adhere to the concept `VolumetricFunction`
 * \ingroup materials
 */
template <Concepts::VolumetricFunction VF>
struct Volumetric
{
  using ScalarType = typename VF::ScalarType;
  using JType      = typename VF::JType;

  using VolumetricFunction = VF;
  using MaterialParameter  = BulkModulus;

  [[nodiscard]] constexpr static std::string name() noexcept { return "Volumetric function: " + VF::name(); }

  /**
   * \brief Construct a new volumetric part
   *
   * \param K Materialparameter is most often the bulk modulus, but it is also common for modified material laws to use
   * lamés first parameter
   * \param vf the volumetric function
   */
  Volumetric(MaterialParameter K, const VF vf)
      : K_{K},
        volumetricFunction_{vf} {}

  /**
   * \brief Computes stored energy of the volumetric function
   *
   * \param J determinant of the deformation gradient $J = \det\BF$
   * \return ScalarType energy
   */
  ScalarType storedEnergy(const JType& J) const { return K_.K * volumetricFunction_.storedEnergyImpl(J); };

  /**
   * \brief Computes the first derivatives of the energy of the volumetric function w.r.t $J$
   *
   * \param J determinant of the deformation gradient $J = \det\BF$
   * \return ScalarType energy
   */
  ScalarType firstDerivative(const JType& J) const { return K_.K * volumetricFunction_.firstDerivativeImpl(J); }

  /**
   * \brief Computes the second derivatives of the energy of the volumetric function w.r.t $J$
   *
   * \param J determinant of the deformation gradient $J = \det\BF$
   * \return ScalarType energy
   */
  ScalarType secondDerivative(const JType& J) const { return K_.K * volumetricFunction_.secondDerivativeImpl(J); };

  template <typename STO>
  auto rebind() const {
    auto reboundVF = volumetricFunction_.template rebind<STO>();
    return Volumetric<decltype(reboundVF)>(K_, reboundVF);
  }

private:
  MaterialParameter K_;
  VF volumetricFunction_;
};
} // namespace Ikarus::Materials
