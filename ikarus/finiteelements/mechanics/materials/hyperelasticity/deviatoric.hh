// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file deviatoric.hh
 * \brief Implementation of the Deviatoric material model.
 * \ingroup  materials
 */

#pragma once

#include <dune/common/float_cmp.hh>

#include <ikarus/utils/tensorutils.hh>


namespace Ikarus {


template <typename DF>
struct DeviatoricPart
{
  using ScalarType         = typename DF::ScalarType;
  using PrincipalStretches = typename DF::PrincipalStretches;
  using MaterialParameters = typename DF::MaterialParameters;

  using FirstDerivative = typename DF::FirstDerivative;

  static constexpr int worldDimension = 3;

  using StressMatrix   = Eigen::Vector<ScalarType, worldDimension>;
  using MaterialTensor = Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>>;

  [[nodiscard]] constexpr static std::string name() noexcept { return "Deviatoric function: " + DF::name(); }

  DeviatoricPart(const DF df)
      : deviatoricFunction_{df} {}

  ScalarType storedEnergyImpl(const PrincipalStretches& lambdas) const {
    return deviatoricFunction_.storedEnergyImpl(lambdas);
  };

  StressMatrix stressesImpl(const PrincipalStretches& lambda) const {
    auto dWdLambda = deviatoricFunction_.firstDerivativeImpl(lambda);

    // Compute the principal PK2 stresses by dividing by the stretches
    StressMatrix S;
    for (auto k : dimensionRange())
      S[k] = dWdLambda[k] / lambda[k];

    return S;
  }

  MaterialTensor tangentModuliImpl(const PrincipalStretches& lambda) const {
    auto S  = stressesImpl(lambda);
    auto dS = deviatoricFunction_.secondDerivativeImpl(lambda);

    // Konvektive coordinates
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

  template <typename STO>
  auto rebind() const {
    auto reboundDF = deviatoricFunction_.template rebind<STO>();
    return DeviatoricPart<decltype(reboundDF)>{reboundDF};
  }

private:
  DF deviatoricFunction_;

  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(worldDimension); }
};
} // namespace Ikarus