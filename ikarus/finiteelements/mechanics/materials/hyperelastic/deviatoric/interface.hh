// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file interface.hh
 * \brief Implementation of the interface for the deviatoric part of a hyperelastic material.
 * \ingroup  materials
 */

#pragma once

#include <dune/common/float_cmp.hh>

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/concepts.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief This is the interface implementation for the deviatoric part of a hyperelastic material.
 *    It is intended to be used with the hyperelastic material model.
 * \ingroup materials
 *
 * \details The deviatoric part of the hyperelastic model, i.e., related to
 * \f$ \hat{\Psi}(\lambda_1, \lambda_2, \lambda_3) \f$, is parametrized with a certain
 * deviatoric function (DF) implemented in terms of principal stretches.
 * The three interface functions (energy, streses and tangentModulus) are called with the argument being the principal
 * stretches (\f$ \lambda_i \f$). The underlying deviatoric function must only implement the energy
 * \f$ \hat{\Psi}(\lambda_1, \lambda_2, \lambda_3) \f$ and its first and second derivatives
 * w.r.t the total principal stretches.
 *
 * \tparam DF Deviatoric function.
 */
template <Concepts::DeviatoricFunction DF>
struct Deviatoric
{
  using ScalarType         = typename DF::ScalarType;
  using PrincipalStretches = typename DF::PrincipalStretches;
  using MaterialParameters = typename DF::MaterialParameters;

  using DeviatoricFunction = DF;

  using FirstDerivative  = typename DF::FirstDerivative;
  using SecondDerivative = typename DF::SecondDerivative;

  static constexpr int dim = 3;

  using StressMatrix   = Eigen::Vector<ScalarType, dim>;
  using MaterialTensor = Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>>;

  [[nodiscard]] constexpr static std::string name() noexcept { return "Deviatoric function: " + DF::name(); }

  Deviatoric(const DF df)
      : deviatoricFunction_{df} {}

  /**
   * \brief Returns the material parameters stored in the deviatoric part of the material.
   */
  const MaterialParameters& materialParametersImpl() const { return deviatoricFunction_.materialParametersImpl(); }

  /**
   * \brief Returns the stored energy obtained from the deviatoric function.
   *
   * \param lambdas the principal stretches
   * \return ScalarType the energy
   */
  ScalarType storedEnergy(const PrincipalStretches& lambda) const {
    checkDuplicates(lambda);
    return deviatoricFunction_.storedEnergyImpl(lambda);
  };

  /**
   * \brief Returns the principal PK2 stresses obtained from the first derivative of the deviatoric function.
   *
   * \param lambda the principal stretches
   * \return StressMatrix the stresses in Cartesian coordinate system.
   */
  StressMatrix stresses(const PrincipalStretches& lambda) const {
    checkDuplicates(lambda);
    auto dWdLambda = deviatoricFunction_.firstDerivativeImpl(lambda);

    // Compute the principal PK2 stresses by dividing by the stretches
    StressMatrix S;
    for (auto k : dimensionRange())
      S[k] = dWdLambda[k] / lambda[k];

    return S;
  }

  /**
   * \brief Returns the material tangent modulus obtained from the second derivative of the deviatoric function.
   *
   * \param lambda the principal stretches
   * \return MaterialTensor the tangentModuli in Cartesian coordinate system.
   */
  MaterialTensor tangentModuli(const PrincipalStretches& lambda) const {
    checkDuplicates(lambda);
    auto S  = stresses(lambda);
    auto dS = deviatoricFunction_.secondDerivativeImpl(lambda);

    auto L = MaterialTensor{};
    L.setZero();

    for (auto i : dimensionRange())
      for (auto k : dimensionRange())
        L(i, i, k, k) = 1.0 / (lambda(i) * lambda(k)) * dS(i, k);

    for (auto i : dimensionRange())
      for (auto k : dimensionRange())
        if (i != k) {
          if (Dune::FloatCmp::eq(lambda(i), lambda(k), 1e-8))
            L(i, k, i, k) = 0.5 * (L(i, i, i, i) - L(i, i, k, k));
          else
            L(i, k, i, k) += (S(i) - S(k)) / (pow(lambda(i), 2) - pow(lambda(k), 2));
        }

    return L;
  };

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return The rebound deviatoric part.
   */
  template <typename STO>
  auto rebind() const {
    auto reboundDF = deviatoricFunction_.template rebind<STO>();
    return Deviatoric<decltype(reboundDF)>{reboundDF};
  }

private:
  DF deviatoricFunction_;

  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(dim); }

  /**
   * \brief A function to check if duplicate principal stretches exists.
   *
   * \details The computation of tangentModuli includes an additional term in the denominator (\lambda_b - \lambda_a).
   * This results in nan if \lambda_b = \lambda_a. In the explicit implementation, this is circumvented using the
   * L'Hôpital's rule. However, AutoDiff doesn't see such numerical issues and hence results in nan. The tolerance
   * can be greater than 1e-4, but that leads to inaccurate results.
   *
   * \param lambda Principal stretches.
   * \param tol Tolerance used during the comparison of the principal stretches.
   */
  void checkDuplicates(const PrincipalStretches& lambda, double tol = 1e-4) const {
    if constexpr (not Concepts::AutodiffScalar<ScalarType>)
      return;
    PrincipalStretches sortedLambda = lambda;
    std::ranges::sort(sortedLambda);
    const bool hasDuplicates =
        std::adjacent_find(sortedLambda.begin(), sortedLambda.end(), [&](ScalarType a, ScalarType b) {
          return Dune::FloatCmp::eq(a, b, tol);
        }) != sortedLambda.end();
    if (hasDuplicates)
      DUNE_THROW(Dune::InvalidStateException, "AutoDiff doesn't work if there are duplicate principal stretches.");
  }
};
} // namespace Ikarus::Materials