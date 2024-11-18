// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file ogden.hh
 * \brief Implementation of the regularized Gent material model.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/deviatoricinvariants.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

struct GentMatParameters
{
  double mu;
  double Jm;
};
} // namespace Ikarus

namespace Ikarus::Materials {

/**
 * \brief Implementation of the Gent material model.
 *
 * \tparam ST The scalar type for the strains and stresses,....
 * \tparam n number of ogden parameters
 * \tparam tag type of principal stretch quantity, either total stretches or deviatoric stretches
 * \ingroup materials
 */
template <typename ST>
struct GentT
{
  static constexpr int dim = 3;
  using ScalarType         = ST;
  using PrincipalStretches = Eigen::Vector<ScalarType, dim>;
  using Invariants         = PrincipalStretches;
  using FirstDerivative    = Eigen::Vector<ScalarType, dim>;
  using SecondDerivative   = Eigen::Matrix<ScalarType, dim, dim>;
  using MaterialParameters = GentMatParameters;

  [[nodiscard]] constexpr static std::string name() noexcept { return "Gent"; }

  /**
   * \brief Constructor for GentT
   *
   * \param C material constant
   * \param lambdaM maximum stretch at which the polymer chain locks
   */
  explicit GentT(const MaterialParameters& matPar)
      : matPar_{matPar} {}

  /**
   * \brief Computes the stored energy in the Gent material model.
   * \details Using only the first five terms of the inverse Langevin function.
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  ScalarType storedEnergyImpl(const PrincipalStretches& lambda) const {
    const Invariants& invariants = Impl::invariants(lambda);
    const auto& devInvariants    = DeviatoricInvariants<PrincipalStretches>(lambda);
    const auto W1                = devInvariants.value().first;
    return -(matPar_.mu / 2.0) * matPar_.Jm * log(1.0 - ((W1 - 3.0) / matPar_.Jm));
  }

  /**
   * \brief Computes the first derivative of the stored energy function w.r.t. the total principal stretches
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  FirstDerivative firstDerivativeImpl(const PrincipalStretches& lambda) const {
    const Invariants& invariants = Impl::invariants(lambda);
    auto dWdLambda               = FirstDerivative::Zero().eval();

    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches>(lambda);
    const auto W1             = devInvariants.value().first;
    const auto& dW1dLambda    = devInvariants.firstDerivative().first;
    const auto mu             = matPar_.mu;
    const auto Jm             = matPar_.Jm;

    for (auto k : dimensionRange())
      dWdLambda[k] += (mu * dW1dLambda[k] * Jm) / (2.0 * (Jm - W1) + 6.0);

    return dWdLambda;
  }

  /**
   * \brief Computes the second derivatives of the stored energy function w.r.t. the total principal stretches
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  SecondDerivative secondDerivativeImpl(const PrincipalStretches& lambda) const {
    const Invariants& invariants = Impl::invariants(lambda);
    auto dS                      = SecondDerivative::Zero().eval();

    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches>(lambda);
    const auto W1             = devInvariants.value().first;
    const auto& dW1dLambda    = devInvariants.firstDerivative().first;
    const auto& ddW1dLambda   = devInvariants.secondDerivative().first;
    const auto mu             = matPar_.mu;
    const auto Jm             = matPar_.Jm;

    for (auto i : dimensionRange())
      for (auto j : dimensionRange()) {
        auto factor1 = 1.0 - ((W1 - 3.0) / Jm);
        dS(i, j) += (mu / (2.0 * (lambda[i] * factor1))) *
                    (ddW1dLambda(i, j) + (dW1dLambda[i] * dW1dLambda[j] / (factor1 * Jm)));
        if (i == j)
          dS(i, j) -= (mu / (2.0 * pow(lambda[i], 2) * factor1)) * dW1dLambda[i];
      }

    return dS;
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return GentT<ScalarTypeOther> The rebound Gent material.
   */
  template <typename STO>
  auto rebind() const {
    return GentT<STO>(matPar_);
  }

private:
  MaterialParameters matPar_;

  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(dim); }
};

/**
 * \brief Alias for GentT with double as the default scalar type.
 */
using Gent = GentT<double>;

} // namespace Ikarus::Materials
