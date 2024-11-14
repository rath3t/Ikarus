// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file BlatzKo.hh
 * \brief Implementation of the BlatzKo material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Implementation of the Blatz-Ko material model.
 *
 * \tparam ST The scalar type for the strains and stresses,....
 */
template <typename ST>
struct BlatzKoT
{
  using ScalarType         = ST;
  using PrincipalStretches = Eigen::Vector<ScalarType, 3>;

  static constexpr int dim         = 3;
  static constexpr auto stretchTag = PrincipalStretchTag::total;

  using FirstDerivative  = Eigen::Vector<ScalarType, dim>;
  using SecondDerivative = Eigen::Matrix<ScalarType, dim, dim>;

  using MaterialParameters = ShearModulus;

  [[nodiscard]] constexpr static std::string name() noexcept { return "BlatzKo"; }

  /**
   * \brief Constructor for BlatzKoT.
   * \param mpt material parameters, here the shear modulus mu.
   */
  explicit BlatzKoT(const MaterialParameters& mpt)
      : materialParameter_{mpt} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameter_; }

  /**
   * \brief Computes the stored energy in the BlatzKo material model.
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  ScalarType storedEnergyImpl(const PrincipalStretches& lambda) const {
    return materialParameter_.mu / 2 *
           (1 / pow(lambda[0], 2) + 1 / pow(lambda[1], 2) + 1 / pow(lambda[2], 2) +
            2 * lambda[0] * lambda[1] * lambda[2] - 5);
  }

  /**
   * \brief Computes the first derivative of the stored energy function w.r.t. the total principal stretches
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  FirstDerivative firstDerivativeImpl(const PrincipalStretches& lambda) const {
    auto dWdLambda = FirstDerivative::Zero().eval();

    dWdLambda[0] = materialParameter_.mu * (-2 / pow(lambda[0], 3) + 2 * lambda[1] * lambda[2]) / 2;
    dWdLambda[1] = materialParameter_.mu * (-2 / pow(lambda[1], 3) + 2 * lambda[0] * lambda[2]) / 2;
    dWdLambda[2] = materialParameter_.mu * (-2 / pow(lambda[2], 3) + 2 * lambda[0] * lambda[1]) / 2;

    return dWdLambda;
  }

  /**
   * \brief Computes the second derivatives of the stored energy function w.r.t. the total principal stretches
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  SecondDerivative secondDerivativeImpl(const PrincipalStretches& lambda) const {
    auto ddWdLambda = SecondDerivative::Zero().eval();
    auto mu         = materialParameter_.mu;

    ddWdLambda(0, 0) = -mu * (-2.0 / pow(lambda(0), 3) + 2.0 * lambda(1) * lambda(2)) / (2.0 * pow(lambda(0), 2)) +
                       3.0 * mu / pow(lambda(0), 5);
    ddWdLambda(0, 1) = mu * lambda(2) / lambda(0);
    ddWdLambda(0, 2) = mu * lambda(1) / lambda(0);
    ddWdLambda(1, 0) = mu * lambda(2) / lambda(1);
    ddWdLambda(1, 1) = -mu * (-2.0 / pow(lambda(1), 3) + 2.0 * lambda(0) * lambda(2)) / (2.0 * pow(lambda(1), 2)) +
                       3.0 * mu / pow(lambda(1), 5);
    ddWdLambda(1, 2) = mu * lambda(0) / lambda(1);
    ddWdLambda(2, 0) = mu * lambda(1) / lambda(2);
    ddWdLambda(2, 1) = mu * lambda(0) / lambda(2);
    ddWdLambda(2, 2) = -mu * (-2.0 / pow(lambda(2), 3) + 2.0 * lambda(0) * lambda(1)) / (2.0 * pow(lambda(2), 2)) +
                       3.0 * mu / pow(lambda(2), 5);

    return ddWdLambda;
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

} // namespace Ikarus::Materials
