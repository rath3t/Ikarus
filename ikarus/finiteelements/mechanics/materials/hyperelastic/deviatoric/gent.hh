// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file gent.hh
 * \brief Implementation of the Gent material model.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/deviatoricinvariants.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

///< Structure representing material parameters for the Gent material model.
struct GentMatParameters
{
  double mu; ///< Denotes the shear modulus.
  double Jm; ///< Denotes a dimensionless parameter that controls the extensibility of chains
};
} // namespace Ikarus

namespace Ikarus::Materials {

/**
 * \brief Implementation of the Gent material model.
 * \ingroup materials
 *
 * \details The energy is computed as
 * \f[ \hat{\Psi}(\lambda_1, \lambda_2, \lambda_3) = -\frac{\mu}{2} J_m \ln{(1-\frac{W_1 - 3}{J_m})}, \f]
 * with \f$ J_m > W_1-3\f$.
 *
 * The first derivatives w.r.t the total principal stretches are
 * \f[ \fracpt{\Psi}{\lambda_i} = \frac{\mu J_m}{2(J_m - W_1) + 6}\fracpt{W_1}{\lambda_i}. \f]
 *
 * The second derivatives w.r.t the total principal stretches are
 * \f[ \fracpt{^2 \Psi}{\lambda_i\partial\lambda_j} =
 \frac{\mu}{2\alpha} (\fracpt{^2 W_1}{\lambda_i\partial\lambda_j} +
 \frac{1}{\alpha J_m} \fracpt{W1}{\lambda_i} \fracpt{W1}{\lambda_j}) -
 \delta_{ij} \frac{\mu}{2\alpha\lambda_i} \fracpt{W_1}{\lambda_i}, \f]
 * with \f$ \alpha = 1 - \frac{W_1 - 3}{J_m} \f].
 *
 * \remark See \cite bergstromMechanicsSolidPolymers2015 for details on this material. For information on the deviatoric
 * invariant \f$ W_1 \f$, see \file deviatoricinvariants.hh 
 *
 * \tparam ST The underlying scalar type.
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

  static constexpr auto stretchTag = PrincipalStretchTag::deviatoric;

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
   * \brief Returns the material parameters stored in the material
   */
  const MaterialParameters& materialParametersImpl() const { return matPar_; }

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
    checkJm(W1);
    return -(matPar_.mu / 2.0) * matPar_.Jm * log(1.0 - ((W1 - 3.0) / matPar_.Jm));
  }

  /**
   * \brief Computes the first derivative of the stored energy function w.r.t. the total principal stretches.
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
    checkJm(W1);

    for (auto k : dimensionRange())
      dWdLambda[k] += (mu * dW1dLambda[k] * Jm) / (2.0 * (Jm - W1) + 6.0);

    return dWdLambda;
  }

  /**
   * \brief Computes the second derivatives of the stored energy function w.r.t. the total principal stretches.
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
    checkJm(W1);

    for (auto i : dimensionRange())
      for (auto j : dimensionRange()) {
        auto factor1 = 1.0 - ((W1 - 3.0) / Jm);
        dS(i, j) += (mu / (2.0 * factor1)) * (ddW1dLambda(i, j) + (dW1dLambda[i] * dW1dLambda[j] / (factor1 * Jm)));
        if (i == j)
          dS(i, j) -= (mu / (2.0 * lambda[i] * factor1)) * dW1dLambda[i];
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

  void checkJm(ScalarType W1) const {
    if (Dune::FloatCmp::le(matPar_.Jm, static_cast<double>(W1) - 3.0, 1e-14))
      DUNE_THROW(Dune::InvalidStateException, "The material parameter Jm should be greater than (W1 - 3)");
  }
};

/**
 * \brief Alias for GentT with double as the default scalar type.
 */
using Gent = GentT<double>;

} // namespace Ikarus::Materials
