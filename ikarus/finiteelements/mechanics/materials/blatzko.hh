// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file BlatzKo.hh
 * \brief Implementation of the BlatzKo material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

/**
 * \brief Implementation of the Neo-Hookean material model.
* \ingroup materials
*  The energy is computed as
*  \f[ \psi(\BC) = \frac{\mu}{2} (\tr \BC-3- 2 \log \sqrt{\det \BC}) + \frac{\lambda}{2} (\log \sqrt{\det \BC})^2 ,\f]
* where \f$ \BC \f$ denotes the right Cauchy-Green strain tensor.
*
*  The second Piola-Kirchhoff stresses are computed as
*      \f[ \BS(\BC) =\fracpt{\psi(\BC)}{\BC} = \mu (\BI-\BC^{-1}) + \lambda \log \sqrt{\det \BC}  \BC^{-1},\f]
*
* and the material tangent moduli are computed as
*      \f[ \BBC(\BC) =\fracpt{^2\psi(\BC)}{\BC^2} =  \lambda \BC^{-1} \otimes  \BC^{-1} + 2 (\mu- \lambda \log
\sqrt{\det \BC} ) \CI,\f]
*      where \f$ \CI_{IJKL} =  \frac{1}{2}({(\BC^{-1})}^{IK}{(\BC^{-1})}^{JL}+{(\BC^{-1})}^{IL} {(\BC^{-1})}^{JK}).\f$
*
*  \remark See \cite bonet2008nonlinear, Section 6.4.3 for a discussion of this material
 * \tparam ST The scalar type for the strains and stresses,....
 */
template <typename ST>
struct BlatzKoT : public Material<BlatzKoT<ST>>
{
  using ScalarType                    = ST;
  static constexpr int worldDimension = 3;
  using StrainMatrix                  = Eigen::Matrix<ScalarType, worldDimension, worldDimension>;
  using StressMatrix                  = StrainMatrix;
  using MaterialParameters            = LamesFirstParameterAndShearModulus;

  static constexpr auto strainTag              = StrainTags::rightCauchyGreenTensor;
  static constexpr auto stressTag              = StressTags::PK2;
  static constexpr auto tangentModuliTag       = TangentModuliTags::Material;
  static constexpr bool energyAcceptsVoigt     = false;
  static constexpr bool stressToVoigt          = false;
  static constexpr bool stressAcceptsVoigt     = false;
  static constexpr bool moduliToVoigt          = false;
  static constexpr bool moduliAcceptsVoigt     = false;
  static constexpr double derivativeFactorImpl = 2;

  [[nodiscard]] constexpr static std::string nameImpl() noexcept { return "BlatzKo"; }

  /**
   * \brief Constructor for BlatzKoT.
   * \param mpt The Lame's parameters (first parameter and shear modulus).
   */
  explicit BlatzKoT(const MaterialParameters& mpt)
      : materialParameter_{mpt} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameter_; }

  // TODO: Can we use SelfAdjointSolver?
  template <typename Derived>
  auto principalStretches(const Eigen::MatrixBase<Derived>& C, int options = Eigen::ComputeEigenvectors) const {
    Eigen::SelfAdjointEigenSolver<Derived> eigensolver(C, options);
    auto& eigenvalues  = eigensolver.eigenvalues();
    auto& eigenvectors = options == Eigen::ComputeEigenvectors ? eigensolver.eigenvectors() : Derived::Zero();

    auto principalStretches = eigenvalues.array().sqrt().eval();
    return std::make_pair(principalStretches, eigenvectors);
  }

  /**
   * \brief Computes the stored energy in the Neo-Hookean material model.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return ScalarType The stored energy.
   */
  ScalarType storedEnergyImpl(const Eigen::Vector<ScalarType, 3>& lambdas) const {
    return materialParameter_.mu / 2 *
           (1 / std::pow(lambdas[0], 2) + 1 / std::pow(lambdas[1], 2) + 1 / std::pow(lambdas[2], 2) +
            2 * lambdas[0] * lambdas[1] * lambdas[2] - 5);
  }

  /**
   * \brief Computes the stresses in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return StressMatrix The stresses.
   */
  auto stressesImpl(const Eigen::Vector<ScalarType, 3>& lambdas) const { return principalStresseses(lambdas); }

  /**
   * \brief Computes the tangent moduli in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> The tangent moduli.
   */

  auto tangentModuliImpl(const Eigen::Vector<ScalarType, 3>& lambdas) const {
    auto S            = principalStresseses(lambdas);
    auto dS = dSdLambda(lambdas);

    // Konvektive coordinates
    auto L = Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>>{};
    L.setZero();

    for (int i = 0; i < 3; ++i) {
      for (int k = 0; k < 3; ++k) {
        L(i, i, k, k) = 1.0 / lambdas(k) * dS(i, k);
      }
    }

    for (int i = 0; i < 3; ++i) {
      for (int k = 0; k < 3; ++k) {
        if (i != k) {
          if (Dune::FloatCmp::eq(lambdas(i), lambdas(k), 1e-8)) {
            L(i, k, i, k) = 0.5 * (L(i, i, i, i) - L(i, i, k, k));
          } else {
            L(i, k, i, k) += (S(i) - S(k)) / (std::pow(lambdas(i), 2) - std::pow(lambdas(k), 2));
          }
        }
      }
    }

    return L;
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return BlatzKoT<ScalarTypeOther> The rebound BlatzKo material.
   */
  template <typename STO>
  auto rebind() const {
    return BlatzKoT<STO>(materialParameter_);
  }

private:
  MaterialParameters materialParameter_;

  auto principalStresseses(const auto& lambdas) const {
    // principal PK1 stress
    auto P1 = materialParameter_.mu * (-2 / std::pow(lambdas[0], 3) + 2 * lambdas[1] * lambdas[2]) / 2;
    auto P2 = materialParameter_.mu * (-2 / std::pow(lambdas[1], 3) + 2 * lambdas[0] * lambdas[2]) / 2;
    auto P3 = materialParameter_.mu * (-2 / std::pow(lambdas[2], 3) + 2 * lambdas[0] * lambdas[1]) / 2;

    // principal PK2 stress
    auto S1 = 1 / lambdas[0] * P1;
    auto S2 = 1 / lambdas[1] * P2;
    auto S3 = 1 / lambdas[2] * P3;

    return Eigen::Vector<ScalarType, 3>{S1, S2, S3};
  }

  auto dSdLambda(const auto& lambda) const {
    auto dS = Eigen::Matrix3<ScalarType>::Zero().eval();

    double mu = materialParameter_.mu;
    dS(0, 0)  = -mu * (-2.0 / std::pow(lambda(0), 3) + 2.0 * lambda(1) * lambda(2)) / (2.0 * std::pow(lambda(0), 2)) +
               3.0 * mu / std::pow(lambda(0), 5);
    dS(0, 1) = mu * lambda(2) / lambda(0);
    dS(0, 2) = mu * lambda(1) / lambda(0);
    dS(1, 0) = mu * lambda(2) / lambda(1);
    dS(1, 1) = -mu * (-2.0 / std::pow(lambda(1), 3) + 2.0 * lambda(0) * lambda(2)) / (2.0 * std::pow(lambda(1), 2)) +
               3.0 * mu / std::pow(lambda(1), 5);
    dS(1, 2) = mu * lambda(0) / lambda(1);
    dS(2, 0) = mu * lambda(1) / lambda(2);
    dS(2, 1) = mu * lambda(0) / lambda(2);
    dS(2, 2) = -mu * (-2.0 / std::pow(lambda(2), 3) + 2.0 * lambda(0) * lambda(1)) / (2.0 * std::pow(lambda(2), 2)) +
               3.0 * mu / std::pow(lambda(2), 5);

    return dS;
  }
};

/**
 * \brief Alias for BlatzKoT with double as the default scalar type.
 */
using BlatzKo = BlatzKoT<double>;

} // namespace Ikarus
