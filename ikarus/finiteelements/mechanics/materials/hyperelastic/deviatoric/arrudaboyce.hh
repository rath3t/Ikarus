// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file ogden.hh
 * \brief Implementation of the regularized ArrudaBoyce material model.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/deviatoricinvariants.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

struct ArrudaBoyceMatParameters
{
  double C_;
  double lambdaM_;
};
} // namespace Ikarus

namespace Ikarus::Materials {

/**
 * \brief Implementation of the ArrudaBoyce material model.
 *
 * \tparam ST The scalar type for the strains and stresses,....
 * \tparam n number of ogden parameters
 * \tparam tag type of principal stretch quantity, either total stretches or deviatoric stretches
 * \ingroup materials
 */
template <typename ST>
struct ArrudaBoyceT
{
  static constexpr int dim      = 3;
  using ScalarType              = ST;
  using PrincipalStretches      = Eigen::Vector<ScalarType, dim>;
  using Invariants              = PrincipalStretches;
  using FirstDerivative         = Eigen::Vector<ScalarType, dim>;
  using SecondDerivative        = Eigen::Matrix<ScalarType, dim, dim>;
  using MaterialParameters      = ArrudaBoyceMatParameters;
  static constexpr int numTerms = 5;

  [[nodiscard]] constexpr static std::string name() noexcept { return "ArrudaBoyce"; }

  /**
   * \brief Constructor for ArrudaBoyceT
   *
   * \param C material constant
   * \param lambdaM maximum stretch at which the polymer chain locks
   */
  explicit ArrudaBoyceT(const MaterialParameters& matPar)
      : matPar_{matPar} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  const MaterialParameters& materialParametersImpl() const { return matPar_; }

  /**
   * \brief Computes the stored energy in the ArrudaBoyce material model.
   * \details Using only the first five terms of the inverse Langevin function.
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  ScalarType storedEnergyImpl(const PrincipalStretches& lambda) const {
    const Invariants& invariants = Impl::invariants(lambda);
    ScalarType energy{0.0};
    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches>(lambda);
    auto W1                   = devInvariants.value().first;
    const auto C_             = matPar_.C_;
    const auto lambdaM_       = matPar_.lambdaM_;
    const auto beta           = 1 / pow(lambdaM_, 2.0);

    for (auto i : parameterRange())
      energy += alphas_[i] * pow(beta, i) * (pow(W1, i + 1) - pow(3, i + 1));
    energy *= C_;

    return energy;
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
    auto W1                   = devInvariants.value().first;
    const auto& dW1dLambda    = devInvariants.firstDerivative().first;
    const auto C_             = matPar_.C_;
    const auto lambdaM_       = matPar_.lambdaM_;
    const auto beta           = 1 / pow(lambdaM_, 2.0);

    for (auto j : parameterRange())
      for (auto k : dimensionRange())
        dWdLambda[k] += C_ * alphas_[j] * pow(beta, j) * pow(W1, j) * dW1dLambda[k] * (j + 1);

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
    auto W1                   = devInvariants.value().first;
    const auto& dW1dLambda    = devInvariants.firstDerivative().first;
    const auto& ddW1dLambda   = devInvariants.secondDerivative().first;
    const auto C_             = matPar_.C_;
    const auto lambdaM_       = matPar_.lambdaM_;
    const auto beta           = 1 / pow(lambdaM_, 2.0);

    for (auto p : parameterRange())
      for (auto i : dimensionRange())
        for (auto j : dimensionRange()) {
          auto factor1 = C_ * alphas_[p] * pow(beta, p);
          auto factor2 = pow(W1, p) * ddW1dLambda(i, j) * (p + 1);
          auto factor3 = pow(W1, p - 1) * dW1dLambda[i] * dW1dLambda[j] * p * (p + 1);
          dS(i, j) += factor1 * (factor2 + factor3);
          if (i == j)
            dS(i, j) -= (1.0 / lambda[i]) * factor1 * pow(W1, p) * dW1dLambda[i] * (p + 1);
        }

    return dS;
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return ArrudaBoyceT<ScalarTypeOther> The rebound ArrudaBoyce material.
   */
  template <typename STO>
  auto rebind() const {
    return ArrudaBoyceT<STO>(matPar_);
  }

private:
  MaterialParameters matPar_;
  std::array<double, numTerms> alphas_ = {0.5, 1.0 / 20.0, 11.0 / 1050.0, 19.0 / 7000.0, 519.0 / 673750.0};

  inline auto parameterRange() const { return Dune::Hybrid::integralRange(numTerms); }
  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(dim); }
};

/**
 * \brief Alias for ArrudaBoyceT with double as the default scalar type.
 */
using ArrudaBoyce = ArrudaBoyceT<double>;

} // namespace Ikarus::Materials
