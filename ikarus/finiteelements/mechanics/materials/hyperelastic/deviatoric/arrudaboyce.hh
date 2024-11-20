// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file arrudaboyce.hh
 * \brief Implementation of the ArrudaBoyce material model.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/deviatoricinvariants.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

///< Structure representing material parameters for the Arrudy-Boyce material model.
struct ArrudaBoyceMatParameters
{
  double mu_;      ///< Denotes the shear modulus.
  double lambdaM_; ///< Denotes the maximum (fully extended) stretch that a molecule is exposed to.
};
} // namespace Ikarus

namespace Ikarus::Materials {

/**
 * \brief Implementation of the ArrudaBoyce material model (also referred as Eight-Chain model).
 * \ingroup materials
 *
 * \details The energy is computed as
 * \f[ \hat{\Psi}(\lambda_1, \lambda_2, \lambda_3) = \mu  \sum_{p=0}^4{\alpha_p  \beta^p}  (W_1^{p+1} - 3^{p+1}), \f]
 * with \f$ \beta = \frac{1}{\lambda_m^2} , \alpha_0 = \frac{1}{2}, \alpha_1 = \frac{1}{20}, \alpha_2 = \frac{11}{1050},
 \alpha_3 = \frac{19}{7000}, \alpha_4 = \frac{519}{673750} \f$.
 *
 * The first derivatives w.r.t the total principal stretches are
 * \f[ \fracpt{\Psi}{\lambda_i} = \mu \sum_{p=0}^4{\alpha_p  \beta^p  W_1^p (p+1)
 \fracpt{W_1}{\lambda_i}}. \f]
 *
 * The second derivatives w.r.t the total principal stretches are
 * \f[ \fracpt{^2 \Psi}{\lambda_i\partial\lambda_j} = \mu \sum_{p=0}^4{
  \alpha_p \beta^p (W_1^p (p+1) \fracpt{^2 W_1}{\lambda_i\partial\lambda_j} +
 W_1^{p-1} p (p+1) \fracpt{W_1}{\lambda_i} \fracpt{W_1}{\lambda_j} - \delta_{ij}
 \frac{1}{\lambda_i} W_1^p (p+1) \fracpt{W_1}{\lambda_i}) } \f]
 *
 * \remark See \cite hiermaierStructuresCrashImpact2010 and \cite bergstromMechanicsSolidPolymers2015 for details on
 * this material. For information on the deviatoric invariant \f$ W_1 \f$, see \ref DeviatoricInvariants.
 *
 *\tparam ST The underlying scalar type.
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

  static constexpr auto stretchTag = PrincipalStretchTag::deviatoric;

  [[nodiscard]] constexpr static std::string name() noexcept { return "ArrudaBoyce"; }

  /**
   * \brief Constructor for ArrudaBoyceT
   *
   * \param matPar The material parameters for the Arruda-Boyce material model.
   */
  explicit ArrudaBoyceT(const MaterialParameters& matPar)
      : matPar_{matPar} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return matPar_; }

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
    const auto W1             = devInvariants.value().first;
    const auto mu_            = matPar_.mu_;
    const auto lambdaM_       = matPar_.lambdaM_;
    const auto beta           = 1 / pow(lambdaM_, 2.0);

    for (const auto i : parameterRange())
      energy += alphas_[i] * pow(beta, i) * (pow(W1, i + 1) - pow(3, i + 1));
    energy *= mu_;

    return energy;
  }

  /**
   * \brief Computes the first derivative of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches
   * \return FirstDerivative The first derivatives of the stored energy function w.r.t. the total principal stretches.
   */
  FirstDerivative firstDerivativeImpl(const PrincipalStretches& lambda) const {
    const Invariants& invariants = Impl::invariants(lambda);
    auto dWdLambda               = FirstDerivative::Zero().eval();

    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches>(lambda);
    const auto W1             = devInvariants.value().first;
    const auto& dW1dLambda    = devInvariants.firstDerivative().first;
    const auto mu_            = matPar_.mu_;
    const auto lambdaM_       = matPar_.lambdaM_;
    const auto beta           = 1 / pow(lambdaM_, 2.0);

    for (const auto j : parameterRange())
      for (const auto k : dimensionRange())
        dWdLambda[k] += mu_ * alphas_[j] * pow(beta, j) * pow(W1, j) * dW1dLambda[k] * (j + 1);

    return dWdLambda;
  }

  /**
   * \brief Computes the second derivatives of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches
   * \return SecondDerivative The second derivatives of the stored energy function w.r.t. the total principal stretches.
   */
  SecondDerivative secondDerivativeImpl(const PrincipalStretches& lambda) const {
    const Invariants& invariants = Impl::invariants(lambda);
    auto dS                      = SecondDerivative::Zero().eval();

    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches>(lambda);
    const auto W1             = devInvariants.value().first;
    const auto& dW1dLambda    = devInvariants.firstDerivative().first;
    const auto& ddW1dLambda   = devInvariants.secondDerivative().first;
    const auto mu_            = matPar_.mu_;
    const auto lambdaM_       = matPar_.lambdaM_;
    const auto beta           = 1 / pow(lambdaM_, 2.0);

    for (const auto p : parameterRange())
      for (const auto i : dimensionRange())
        for (const auto j : dimensionRange()) {
          auto factor1 = mu_ * alphas_[p] * pow(beta, p);
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

  inline auto parameterRange() const { return Dune::range(numTerms); }
  inline auto dimensionRange() const { return Dune::range(dim); }
};

/**
 * \brief Alias for ArrudaBoyceT with double as the default scalar type.
 */
using ArrudaBoyce = ArrudaBoyceT<double>;

} // namespace Ikarus::Materials
