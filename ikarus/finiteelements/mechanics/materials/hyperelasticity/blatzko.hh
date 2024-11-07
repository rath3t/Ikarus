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


template <typename ST>
struct BlatzKoT
{
  using ScalarType         = ST;
  using PrincipalStretches = Eigen::Vector<ScalarType, 3>;

  static constexpr int worldDimension = 3;

  using FirstDerivative  = Eigen::Vector<ScalarType, worldDimension>;
  using SecondDerivative = Eigen::Matrix<ScalarType, worldDimension, worldDimension>;

  using MaterialParameters = ShearModulus;

  [[nodiscard]] constexpr static std::string name() noexcept { return "BlatzKo"; }

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


  ScalarType storedEnergyImpl(const PrincipalStretches& lambdas) const {
    return materialParameter_.mu / 2 *
           (1 / pow(lambdas[0], 2) + 1 / pow(lambdas[1], 2) + 1 / pow(lambdas[2], 2) +
            2 * lambdas[0] * lambdas[1] * lambdas[2] - 5);
  }


  FirstDerivative firstDerivativeImpl(const PrincipalStretches& lambdas) const {
    // principal PK1 stress
    auto P1 = materialParameter_.mu * (-2 / pow(lambdas[0], 3) + 2 * lambdas[1] * lambdas[2]) / 2;
    auto P2 = materialParameter_.mu * (-2 / pow(lambdas[1], 3) + 2 * lambdas[0] * lambdas[2]) / 2;
    auto P3 = materialParameter_.mu * (-2 / pow(lambdas[2], 3) + 2 * lambdas[0] * lambdas[1]) / 2;

    return FirstDerivative{P1, P2, P3};
  }

  SecondDerivative secondDerivativeImpl(const PrincipalStretches& lambda) const {
    auto dS = SecondDerivative::Zero().eval();

    double mu = materialParameter_.mu;
    dS(0, 0)  = -mu * (-2.0 / pow(lambda(0), 3) + 2.0 * lambda(1) * lambda(2)) / (2.0 * pow(lambda(0), 2)) +
               3.0 * mu / pow(lambda(0), 5);
    dS(0, 1) = mu * lambda(2) / lambda(0);
    dS(0, 2) = mu * lambda(1) / lambda(0);
    dS(1, 0) = mu * lambda(2) / lambda(1);
    dS(1, 1) = -mu * (-2.0 / pow(lambda(1), 3) + 2.0 * lambda(0) * lambda(2)) / (2.0 * pow(lambda(1), 2)) +
               3.0 * mu / pow(lambda(1), 5);
    dS(1, 2) = mu * lambda(0) / lambda(1);
    dS(2, 0) = mu * lambda(1) / lambda(2);
    dS(2, 1) = mu * lambda(0) / lambda(2);
    dS(2, 2) = -mu * (-2.0 / pow(lambda(2), 3) + 2.0 * lambda(0) * lambda(1)) / (2.0 * pow(lambda(2), 2)) +
               3.0 * mu / pow(lambda(2), 5);

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
