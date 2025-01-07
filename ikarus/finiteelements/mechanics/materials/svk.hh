// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file svk.hh
 * \brief Implementation of the Saint Venant-Kirchhoff material model.
 * \ingroup  materials
 */

// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Implementation of the Saint Venant-Kirchhoff material model.
 * \ingroup materials
 *   The energy is computed as
 *  \f[ \psi(\BE) = \frac{\la}{2} (\tr \BE)^2   +\mu \tr (\BE^2) ,\f]
 *  where \f$ \BE \f$ denotes the Green-Lagrangian strain.
 *
 *  The second Piola-Kirchhoff stresses are computed as
 *     \f[ \BS(\BE) =\fracpt{\psi(\BE)}{\BE} = \la \tr \BE \BI  +2 \mu \BE,\f]
 *
 * and the material tangent moduli are computed as
 *      \f[ \BBC(\BE) =\fracpt{^2\psi(\BE)}{\BE^2} =  \la \tr \BE \CI  +2 \mu \CI^{\mathrm{sym}},\f]
 *      where \f$ \CI_{IJKL} =  \de_{IJ}\de_{KL}\f$ and \f$ \CI_{IJKL}^\mathrm{sym} =  \frac{1}{2}(\de_{IK}\de_{JL}+
 * \de_{IL}\de_{JK})\f$.
 * \tparam ST The scalar type used in the material.
 */
template <typename ST>
struct StVenantKirchhoffT : public Material<StVenantKirchhoffT<ST>>
{
  using ScalarType         = ST;
  static constexpr int dim = 3;
  using StrainMatrix       = Eigen::Matrix<ScalarType, dim, dim>;
  using StressMatrix       = StrainMatrix;
  using MaterialTensor     = Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>>;

  using MaterialParameters = LamesFirstParameterAndShearModulus;

  static constexpr auto strainTag              = StrainTags::greenLagrangian;
  static constexpr auto stressTag              = StressTags::PK2;
  static constexpr auto tangentModuliTag       = TangentModuliTags::Material;
  static constexpr bool energyAcceptsVoigt     = true;
  static constexpr bool stressToVoigt          = true;
  static constexpr bool stressAcceptsVoigt     = true;
  static constexpr bool moduliToVoigt          = true;
  static constexpr bool moduliAcceptsVoigt     = true;
  static constexpr double derivativeFactorImpl = 1;

  [[nodiscard]] constexpr static std::string nameImpl() { return "StVenantKirchhoff"; }

  /**
   * \brief Constructor for StVenantKirchhoffT.
   * \param mpt The material parameters (Lamé's first parameter and shear modulus).
   */
  explicit StVenantKirchhoffT(const MaterialParameters& mpt)
      : materialParameter_{mpt} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameter_; }

  /**
   * \brief Computes the stored energy in the Saint Venant-Kirchhoff material model.
   * \tparam Derived The derived type of the input matrix.
   * \param E The Green-Lagrangian strain.
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& E) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (Concepts::EigenVector<Derived>) {
      const ScalarType traceE = E.template segment<3>(0).sum();
      const ScalarType squaredNorm =
          E.template segment<3>(0).squaredNorm() + E.template segment<3>(3).squaredNorm() / ScalarType(2.0);
      return materialParameter_.lambda / ScalarType(2.0) * traceE * traceE + materialParameter_.mu * squaredNorm;
    } else {
      const auto traceE = E.trace();
      return materialParameter_.lambda / ScalarType(2.0) * traceE * traceE + materialParameter_.mu * E.squaredNorm();
    }
  }

  /**
   * \brief Computes the stresses in the Saint Venant-Kirchhoff material model.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The Green-Lagrangian strain.
   * \return StressMatrix The stresses.
   */
  template <bool voigt, typename Derived>
  auto stressesImpl(const Eigen::MatrixBase<Derived>& E) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<decltype(E.eval())>);
    const auto& Ed = E.derived();
    if constexpr (!voigt) {
      if constexpr (Concepts::EigenVector<Derived>) {
        static_assert(Concepts::EigenVector6<Derived>);
        StressMatrix S;
        const ScalarType traceE = Ed.template segment<3>(0).sum();
        S.diagonal().array() =
            materialParameter_.lambda * traceE + 2 * materialParameter_.mu * Ed.template segment<3>(0).array();
        S(1, 2) = S(2, 1) = materialParameter_.mu * Ed(3);
        S(0, 2) = S(2, 0) = materialParameter_.mu * Ed(4);
        S(0, 1) = S(1, 0) = materialParameter_.mu * Ed(5); // no two since E voigt has 2* on off-diagonal terms
        return S;
      } else {
        static_assert(Concepts::EigenMatrix33<Derived>);
        return (materialParameter_.lambda * Ed.trace() * StrainMatrix::Identity() + 2 * materialParameter_.mu * Ed)
            .eval();
      }
    } else {
      if constexpr (Concepts::EigenVector<Derived>) {
        static_assert(Concepts::EigenVector6<Derived>);
        Eigen::Matrix<ScalarType, 6, 1> S;
        const ScalarType traceE          = Ed.template segment<3>(0).sum();
        S.template segment<3>(0).array() = traceE * materialParameter_.lambda;
        S.template segment<3>(0) += materialParameter_.mu * 2 * Ed.template segment<3>(0);
        S.template segment<3>(3) =
            materialParameter_.mu * Ed.template segment<3>(3); // no two since E voigt has 2* on off-diagonal terms
        return S;
      } else {
        Eigen::Matrix<ScalarType, 6, 1> S;
        S.template segment<3>(0).array() = Ed.trace() * materialParameter_.lambda;
        S.template segment<3>(0) += 2 * materialParameter_.mu * Ed.diagonal();
        S(3) = 2 * materialParameter_.mu * Ed(1, 2);
        S(4) = 2 * materialParameter_.mu * Ed(0, 2);
        S(5) = 2 * materialParameter_.mu * Ed(0, 1);
        return S;
      }
    }
  }

  /**
   * \brief Computes the tangent moduli in the Saint Venant-Kirchhoff material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The Green-Lagrangian strain (not used).
   * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> or Eigen::Matrix<ScalarType, 6, 6> The tangent
   * moduli.
   */
  template <bool voigt, typename Derived>
  auto tangentModuliImpl([[maybe_unused]] const Eigen::MatrixBase<Derived>& /* E */) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      MaterialTensor moduli;
      moduli = materialParameter_.lambda * identityFourthOrder() +
               2 * materialParameter_.mu * symmetricIdentityFourthOrder();
      return moduli;
    } else {
      Eigen::Matrix<ScalarType, 6, 6> moduli;
      moduli.setZero();
      moduli.template block<3, 3>(0, 0).array() = materialParameter_.lambda;
      moduli.template block<3, 3>(0, 0).diagonal().array() += 2 * materialParameter_.mu;
      moduli.template block<3, 3>(3, 3).diagonal().array() = materialParameter_.mu;
      return moduli;
    }
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam ScalarTypeOther The target scalar type.
   * \return StVenantKirchhoffT<ScalarTypeOther> The rebound StVenantKirchhoff material.
   */
  template <typename ScalarTypeOther>
  auto rebind() const {
    return StVenantKirchhoffT<ScalarTypeOther>(materialParameter_);
  }

private:
  MaterialParameters materialParameter_;
};

/**
 * \brief Alias for StVenantKirchhoffT with double as the default scalar type.
 */
using StVenantKirchhoff = StVenantKirchhoffT<double>;

} // namespace Ikarus::Materials
