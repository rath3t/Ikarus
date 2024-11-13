// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file ogden.hh
 * \brief Implementation of the regularized InvariantBased material model.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Implementation of the InvariantBased material model.
 *
 * \tparam ST The scalar type for the strains and stresses,....
 * \tparam n number of ogden parameters
 * \tparam tag type of principal stretch quantity, either total stretches or deviatoric stretches
 * \ingroup materials
 */
template <typename ST, int n>
struct InvariantBasedT
{
  static constexpr int dim              = 3;
  using ScalarType                      = ST;
  using PrincipalStretches              = Eigen::Vector<ScalarType, dim>;
  using Invariants                      = PrincipalStretches;
  using InvariantReduced                = Eigen::Vector<ScalarType, dim - 1>;
  static constexpr int numMatParameters = n;

  using ExponentsP         = std::array<std::size_t, numMatParameters>;
  using ExponentsQ         = std::array<std::size_t, numMatParameters>;
  using MaterialParameters = std::array<double, numMatParameters>;

  using FirstDerivative  = Eigen::Vector<ScalarType, dim>;
  using SecondDerivative = Eigen::Matrix<ScalarType, dim, dim>;

  [[nodiscard]] constexpr static std::string name() noexcept { return "InvariantBased"; }

  /**
   * \brief Constructor for InvariantBasedT
   *
   * \param mpt material parameters (array of mu values)
   * \param opt ogden parameters (array of alpha values)
   */
  explicit InvariantBasedT(const ExponentsP& pex, const ExponentsQ& qex, const MaterialParameters& mat)
      : pex_{pex},
        qex_{qex},
        mat_{mat} {}

  /**
   * \brief Computes the stored energy in the InvariantBased material model.
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  ScalarType storedEnergyImpl(const PrincipalStretches& lambda) const {
    InvariantReduced I = DeviatoricInvariants(lambda);
    ScalarType energy{};

    for (auto i : parameterRange())
      energy += mat_[i] * pow(I[0], pex_[i]) * pow(I[1], qex_[i]);

    return energy;
  }

  /**
   * \brief Computes the first derivative of the stored energy function w.r.t. the total principal stretches
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  FirstDerivative firstDerivativeImpl(const PrincipalStretches& lambda) const {
    InvariantReduced I = DeviatoricInvariants(lambda);
    auto dWdLambda     = FirstDerivative::Zero().eval();

    return dWdLambda;
  }

  /**
   * \brief Computes the second derivatives of the stored energy function w.r.t. the total principal stretches
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  SecondDerivative secondDerivativeImpl(const PrincipalStretches& lambda) const {
    InvariantReduced I = DeviatoricInvariants(lambda);
    auto dS            = SecondDerivative::Zero().eval();

    return dS;
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return InvariantBasedT<ScalarTypeOther> The rebound InvariantBased material.
   */
  template <typename STO>
  auto rebind() const {
    return InvariantBasedT<STO, numMatParameters>(pex_, qex_, mat_);
  }

private:
  ExponentsP pex_;
  ExponentsQ qex_;
  MaterialParameters mat_;

  InvariantReduced DeviatoricInvariants(const PrincipalStretches& lambda) {
    Invariants invariants = Impl::invariants(lambda);
    auto I                = InvariantReduced::Zero().eval();
    I[0]                  = invariants[0] * pow(invariants[2], -1.0 / 3.0) - 3.0;
    I[1]                  = invariants[1] * pow(invariants[2], -2.0 / 3.0) - 3.0;
    return I;
  }

  inline auto parameterRange() const { return Dune::Hybrid::integralRange(numMatParameters); }
};

/**
 * \brief Alias for InvariantBasedT with double as the default scalar type.
 */
template <int n>
using InvariantBased = InvariantBasedT<double, n>;

} // namespace Ikarus::Materials
