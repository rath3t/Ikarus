// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file deviatoric.hh
 * \brief Implementation of the Deviatoric material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

template <typename DF, bool useIsochoricStretches = true>
struct DeviatoricPart
{
  using ScalarType         = typename DF::ScalarType;
  using PrincipalStretches = typename DF::PrincipalStretches;
  using MaterialParameters = typename DF::MaterialParameters;

  static constexpr int worldDimension = 3;

  using StressMatrix   = Eigen::Vector<ScalarType, worldDimension>;
  using MaterialTensor = Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>>;

  DeviatoricPart(const DF df)
      : deviatoricFunction_{df} {}

  ScalarType storedEnergyImpl(const PrincipalStretches& lambdas) const {
    return deviatoricFunction_.storedEnergyImpl(transformStretches(lambdas));
  };

  StressMatrix stressesImpl(const PrincipalStretches& lambda) const {
    auto lambdaBar = transformStretches(lambda);
    auto P = deviatoricFunction_.firstDerivativeImpl(lambdaBar);

    // Compute the principal PK2 stresses by dividing by the stretches

    StressMatrix S;
    for (auto k : dimensionRange()) {
      S[k] = P[k] / lambdaBar[k];
    }

    return S;
  }

  MaterialTensor tangentModuliImpl(const PrincipalStretches& lambda) const {
    auto lambdaBar = transformStretches(lambda);

    auto S  = stressesImpl(lambdaBar);
    auto dS = deviatoricFunction_.secondDerivativeImpl(lambdaBar);

    // Konvektive coordinates
    auto L = MaterialTensor{};
    L.setZero();

    for (auto i : dimensionRange()) {
      for (auto k : dimensionRange()) {
        L(i, i, k, k) = 1.0 / lambdaBar(k) * dS(i, k);
      }
    }

    for (auto i : dimensionRange()) {
      for (auto k : dimensionRange()) {
        if (i != k) {
          if (Dune::FloatCmp::eq(lambdaBar(i), lambdaBar(k), 1e-8)) {
            L(i, k, i, k) = 0.5 * (L(i, i, i, i) - L(i, i, k, k));
          } else {
            L(i, k, i, k) += (S(i) - S(k)) / (std::pow(lambdaBar(i), 2) - std::pow(lambdaBar(k), 2));
          }
        }
      }
    }

    return L;
  };

private:
  DF deviatoricFunction_;

  inline auto dimensionRange() const { return Dune::Hybrid::integralRange(worldDimension); }

  PrincipalStretches transformStretches(const PrincipalStretches& lambdas) const {
    if constexpr (useIsochoricStretches) {
      ScalarType J    = std::accumulate(lambdas.begin(), lambdas.end(), 1.0, std::multiplies());
      ScalarType Jmod = std::pow(J, -1.0 / 3.0);

      auto lambdasBar = PrincipalStretches::Zero().eval();
      for (auto i : dimensionRange())
        lambdasBar[i] = Jmod * lambdas[i];

      return lambdasBar;
    } else
      return lambdas;
  }
};
} // namespace Ikarus