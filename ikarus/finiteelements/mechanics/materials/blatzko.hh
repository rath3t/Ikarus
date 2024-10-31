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
  using ScalarType         = ST;
  using PrincipalStretches = Eigen::Vector<ScalarType, 3>;

  static constexpr int worldDimension = 3;

  using FirstDerivative  = Eigen::Vector<ScalarType, worldDimension>;
  using SecondDerivative = Eigen::Matrix<ScalarType, worldDimension, worldDimension>;

  using MaterialParameters = ShearModulus;

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

  /**
   * \brief Computes the stored energy in the Neo-Hookean material model.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return ScalarType The stored energy.
   */
  ScalarType storedEnergyImpl(const PrincipalStretches& lambdas) const {
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
  FirstDerivative firstDerivativeImpl(const PrincipalStretches& lambdas) const {
    // principal PK1 stress
    auto P1 = materialParameter_.mu * (-2 / std::pow(lambdas[0], 3) + 2 * lambdas[1] * lambdas[2]) / 2;
    auto P2 = materialParameter_.mu * (-2 / std::pow(lambdas[1], 3) + 2 * lambdas[0] * lambdas[2]) / 2;
    auto P3 = materialParameter_.mu * (-2 / std::pow(lambdas[2], 3) + 2 * lambdas[0] * lambdas[1]) / 2;

    return FirstDerivative{P1, P2, P3};
  }

  /**
   * \brief Computes the tangent moduli in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> The tangent moduli.
   */

  SecondDerivative secondDerivativeImpl(const PrincipalStretches& lambda) const {
    auto dS = SecondDerivative::Zero().eval();

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
  
};

/**
 * \brief Alias for BlatzKoT with double as the default scalar type.
 */
using BlatzKo = BlatzKoT<double>;

} // namespace Ikarus
