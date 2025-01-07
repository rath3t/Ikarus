// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file linearelasticity.hh
 * \brief Contains the LinearElasticityT class implementation.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/svk.hh>

namespace Ikarus::Materials {

/**
 * \brief Implementation of the Linear Elasticity material model.
 * \ingroup materials
 *   The energy is computed as
 *  \f[ \psi(\Bvep) = \frac{\la}{2} (\tr \Bvep)^2   +\mu \tr (\Bvep^2) ,\f]
 *  where \f$ \Bvep \f$ denotes the linear strain tensor.
 *
 *  The second Piola-Kirchhoff stresses are computed as
 *     \f[ \BS(\Bvep) =\fracpt{\psi(\Bvep)}{\Bvep} = \la \tr \Bvep \BI  +2 \mu \Bvep,\f]
 *
 * and the material tangent moduli are computed as
 *      \f[ \BBC(\Bvep) =\fracpt{^2\psi(\Bvep)}{\Bvep^2} =  \la \tr \Bvep \CI  +2 \mu \CI^{\mathrm{sym}},\f]
 *      where \f$ \CI_{IJKL} =  \de_{IJ}\de_{KL}\f$ and \f$ \CI_{IJKL}^\mathrm{sym} =  \frac{1}{2}(\de_{IK}\de_{JL}+
 * \de_{IL}\de_{JK})\f$.
 * \tparam ST The scalar type used in the material.
 */
template <typename ST>
struct LinearElasticityT : Material<LinearElasticityT<ST>>
{
  using ScalarType = ST;
  using Base       = StVenantKirchhoffT<ScalarType>;

  using field_type         = ScalarType;
  static constexpr int dim = Base::dim;
  using StrainMatrix       = typename Base::StrainMatrix;
  using StressMatrix       = StrainMatrix;

  using MaterialParameters = typename Base::MaterialParameters;

  static constexpr auto strainTag              = StrainTags::linear;
  static constexpr auto stressTag              = StressTags::linear;
  static constexpr auto tangentModuliTag       = TangentModuliTags::Material;
  static constexpr bool energyAcceptsVoigt     = Base::energyAcceptsVoigt;
  static constexpr bool stressToVoigt          = Base::stressToVoigt;
  static constexpr bool stressAcceptsVoigt     = Base::stressAcceptsVoigt;
  static constexpr bool moduliToVoigt          = Base::moduliToVoigt;
  static constexpr bool moduliAcceptsVoigt     = Base::moduliAcceptsVoigt;
  static constexpr double derivativeFactorImpl = Base::derivativeFactorImpl;

  [[nodiscard]] constexpr static std::string nameImpl() noexcept { return "LinearElasticity"; }

  /**
   * \brief Constructor for LinearElasticityT.
   *
   * \param mpt The LamesFirstParameterAndShearModulus object.
   */
  explicit LinearElasticityT(const MaterialParameters& mpt)
      : svk_{mpt} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return svk_.materialParametersImpl(); }

  /**
   * \brief Calculates the stored energy in the material.
   *
   * \tparam Derived The underlying Eigen type.
   * \param E The strain tensor components.
   * \return Scalar return of stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& E) const {
    return svk_.template storedEnergyImpl(E);
  }

  /**
   * \brief Calculates the stresses in the material.
   *
   * \tparam voigt Boolean indicating whether to return Voigt-shaped result.
   * \tparam Derived The underlying Eigen type.
   * \param E The strain tensor components.
   * \return Matrix valued return of stresses.
   */
  template <bool voigt, typename Derived>
  auto stressesImpl(const Eigen::MatrixBase<Derived>& E) const {
    return svk_.template stressesImpl<voigt>(E);
  }

  /**
   * \brief Calculates the tangent moduli in the material.
   *
   * \tparam voigt Boolean indicating whether to return Voigt-shaped result.
   * \tparam Derived The underlying Eigen type.
   * \param E The strain tensor components.
   * \return Tangent moduli as fourth-order tensor.
   */
  template <bool voigt, typename Derived>
  auto tangentModuliImpl(const Eigen::MatrixBase<Derived>& E) const {
    return svk_.template tangentModuliImpl<voigt>(E);
  }

  /**
   * \brief Rebind material to a different scalar type.
   *
   * \tparam STO The scalar type to rebind to.
   * \return Rebound material.
   */
  template <typename STO>
  auto rebind() const {
    LinearElasticityT<STO> leRebound;
    leRebound.svk = svk_.template rebind<STO>();
    return leRebound;
  }

private:
  StVenantKirchhoffT<ScalarType> svk_;
};

/**
 * \brief Convenience typedef for LinearElasticity with double as ScalarType.
 */
using LinearElasticity = LinearElasticityT<double>;

} // namespace Ikarus::Materials
