// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file Ogden.hh
 * \brief Implementation of the Ogden material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

template <typename ST, int n>
struct OgdenT : public Material<OgdenT<ST, n>>
{
  using ScalarType                      = ST;
  static constexpr int nOgdenParameters = n;

  static constexpr int worldDimension = 3;

  using PrincipalStretches = Eigen::Vector<ScalarType, 3>;
  using StressMatrix       = Eigen::Vector<ScalarType, worldDimension>;
  using MaterialTensor     = Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>>;

  using MaterialParameters = std::array<double, nOgdenParameters>;
  using OgdenParameters    = std::array<double, nOgdenParameters>;

  static constexpr auto strainTag              = StrainTags::rightCauchyGreenTensor;
  static constexpr auto stressTag              = StressTags::PK2;
  static constexpr auto tangentModuliTag       = TangentModuliTags::Material;
  static constexpr bool energyAcceptsVoigt     = false;
  static constexpr bool stressToVoigt          = false;
  static constexpr bool stressAcceptsVoigt     = false;
  static constexpr bool moduliToVoigt          = false;
  static constexpr bool moduliAcceptsVoigt     = false;
  static constexpr double derivativeFactorImpl = 2;

  [[nodiscard]] constexpr static std::string nameImpl() noexcept { return "Ogden"; }

  /**
   * \brief Constructor for OgdenT.
   * \param mpt The Lame's parameters (first parameter and shear modulus).
   */
  explicit OgdenT(const MaterialParameters& mpt, const OgdenParameters& opt)
      : materialParameters_{mpt},
        ogdenParameters_{opt} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameters_; }

  /**
   * \brief Computes the stored energy in the Neo-Hookean material model.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return ScalarType The stored energy.
   */
  ScalarType storedEnergyImpl(const PrincipalStretches& lambdas) const {
    ScalarType energy{};

    for (auto i : parameterRange()) {
      auto mu_i    = materialParameters_[i];
      auto alpha_i = ogdenParameters_[i];

      energy += mu_i / alpha_i *
                (std::pow(lambdas[0], alpha_i), std::pow(lambdas[1], alpha_i) + std::pow(lambdas[2], alpha_i) - 3);
    }

    return energy;
  }

  /**
   * \brief Computes the stresses in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return StressMatrix The stresses.
   */
  StressMatrix stressesImpl(const PrincipalStretches& lambdas) const { return principalStresseses(lambdas); }

  /**
   * \brief Computes the tangent moduli in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> The tangent moduli.
   */

  MaterialTensor tangentModuliImpl(const PrincipalStretches& lambdas) const {
    auto S  = principalStresseses(lambdas);
    auto dS = dSdLambda(lambdas);

    // Konvektive coordinates
    auto L = MaterialTensor{};
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
   * \return OgdenT<ScalarTypeOther> The rebound Ogden material.
   */
  template <typename STO>
  auto rebind() const {
    return OgdenT<STO, nOgdenParameters>(materialParameters_);
  }

private:
  MaterialParameters materialParameters_;
  OgdenParameters ogdenParameters_;

  inline auto parameterRange() const { return Dune::range(ogdenParameters_.size()); }
  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(worldDimension); }

  auto principalStresseses(const auto& lambdas) const {
    // principal PK1 stress
    auto P = Eigen::Vector3<ScalarType>::Zero().eval();

    auto mu    = materialParameters_;
    auto alpha = ogdenParameters_;

    for (auto j : parameterRange())
      for (auto k : dimensionRange()) {
        P[k] += mu[j] * std::pow(lambdas[k], alpha[j] - 1);
      }

    // principal PK2 stress
    auto S1 = 1 / lambdas[0] * P[0];
    auto S2 = 1 / lambdas[1] * P[1];
    auto S3 = 1 / lambdas[2] * P[2];

    return Eigen::Vector<ScalarType, 3>{S1, S2, S3};
  }

  auto dSdLambda(const auto& lambdas) const {
    auto dS = Eigen::Matrix3<ScalarType>::Zero().eval();

    auto mu    = materialParameters_;
    auto alpha = ogdenParameters_;

    for (auto j : parameterRange())
      for (auto k : dimensionRange()) {
        dS(k, k) += mu[j] * (alpha[j] - 1) * std::pow(lambdas[k], alpha[j] - 2);
      }

    return dS;
  }
};

/**
 * \brief Alias for OgdenT with double as the default scalar type.
 */
template <int n>
using Ogden = OgdenT<double, n>;

} // namespace Ikarus
