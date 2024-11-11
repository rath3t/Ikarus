// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file deviatoric.hh
 * \brief Implementation of the deviatoric part of a hyperelastic material.
 * \ingroup  materials
 */

#pragma once

#include <dune/common/float_cmp.hh>

#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/concepts.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

/**
 * \brief This is the interface implementation for the deviatoric part of a hyperelastic material.
 *    It is intended to use with the hyperelastic material model.
 * \details
 * The deviatoric part is parametrized with a certain deviatoric function (DF) implemented in terms of principal
 * stretches. The three interface functions (energy, streses and tangentModulus) are called with the principal stretches
 * lambda. After calling the deviaoric funtion certain transformation happen to yield the principal stresses and the
 * material tangent in principal coordinates
 *
 * \tparam DF deviatoric material function
 * \ingroup materials
 */
template <Concepts::DeviatoricFunction DF>
struct Deviatoric
{
  using ScalarType         = typename DF::ScalarType;
  using PrincipalStretches = typename DF::PrincipalStretches;
  using MaterialParameters = typename DF::MaterialParameters;

  using DeviatoricFunction = DF;

  using FirstDerivative = typename DF::FirstDerivative;

  static constexpr int dim = 3;

  using StressMatrix   = Eigen::Vector<ScalarType, dim>;
  using MaterialTensor = Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>>;

  [[nodiscard]] constexpr static std::string name() noexcept { return "Deviatoric function: " + DF::name(); }

  Deviatoric(const DF df)
      : deviatoricFunction_{df} {}

  /**
   * \brief Returns the stored energy obtained from the deviatoric function
   *
   * \param lambdas the principal stretches
   * \return ScalarType the energy
   */
  ScalarType storedEnergy(const PrincipalStretches& lambdas) const {
    return deviatoricFunction_.storedEnergyImpl(lambdas);
  };

  /**
   * \brief Returns the principal PK2 stresses obtained from the first derivative of the deviatoric function
   *
   * \param lambda the principal stretches
   * \return StressMatrix
   */
  StressMatrix stresses(const PrincipalStretches& lambda) const {
    auto dWdLambda = deviatoricFunction_.firstDerivativeImpl(lambda);

    // Compute the principal PK2 stresses by dividing by the stretches
    StressMatrix S;
    for (auto k : dimensionRange())
      S[k] = dWdLambda[k] / lambda[k];

    return S;
  }

  /**
   * \brief Returns the material tangent modulus obtained from the second derivative of the deviatoric function
   *
   * \param lambda the principal stretches
   * \return MaterialTensor
   */
  MaterialTensor tangentModuli(const PrincipalStretches& lambda) const {
    auto S  = stresses(lambda);
    auto dS = deviatoricFunction_.secondDerivativeImpl(lambda);

    auto L = MaterialTensor{};
    L.setZero();

    for (auto i : dimensionRange())
      for (auto k : dimensionRange())
        L(i, i, k, k) = 1.0 / lambda(k) * dS(i, k);

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
};
} // namespace Ikarus