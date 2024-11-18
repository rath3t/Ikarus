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
    auto dWdLambda     = FirstDerivative::Zero().eval();
    const ScalarType J = lambda[0] * lambda[1] * lambda[2];

    for (auto k : dimensionRange())
      dWdLambda[k] = materialParameter_.mu * (-2 / pow(lambda[k], 3) + 2 * (J / lambda[k])) / 2;

    return dWdLambda;
  }

  /**
   * \brief Computes the second derivatives of the stored energy function w.r.t. the total principal stretches
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  SecondDerivative secondDerivativeImpl(const PrincipalStretches& lambda) const {
    auto dS            = SecondDerivative::Zero().eval();
    const ScalarType J = lambda[0] * lambda[1] * lambda[2];

    for (auto i : dimensionRange())
      for (auto j : dimensionRange()) {
        if (i == j)
          dS(i, j) += (1.0 / pow(lambda[i], 2)) * (1.0 / pow(lambda[i], 2) - J) + 3.0 / pow(lambda(i), 4);
        else
          dS(i, j) += J / (lambda[i] * lambda[j]);
        dS(i, j) *= materialParameter_.mu;
      }

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

  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(dim); }
};

/**
 * \brief Alias for BlatzKoT with double as the default scalar type.
 */
using BlatzKo = BlatzKoT<double>;

} // namespace Ikarus::Materials
