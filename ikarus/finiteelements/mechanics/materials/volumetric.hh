// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file BlatzKo.hh
 * \brief Implementation of the BlatzKo material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

template <typename VF>
struct VolumetricPart
{
  using ScalarType = typename VF::ScalarType;
  using JType      = typename VF::JType;

  VolumetricPart(BulkModulus K, const VF vf)
      : K_{K},
        volumetricFunction_{vf} {}

  ScalarType storedEnergy(const JType& J) const { return K_.K * volumetricFunction_.storedEnergyImpl(J); };

  ScalarType firstDerivative(const JType& J) const { return K_.K * volumetricFunction_.firstDerivativeImpl(J); }

  ScalarType secondDerivative(const JType& J) const { return K_.K * volumetricFunction_.secondDerivativeImpl(J); };

private:
  BulkModulus K_;
  VF volumetricFunction_;
};

template <typename ST>
struct VF1T
{
  using ScalarType = ST;
  using JType      = ScalarType;

  ScalarType storedEnergyImpl(const JType& J) const { return 0.5 * std::pow(J - 1, 2); };

  ScalarType firstDerivativeImpl(const JType& J) const { return J - 1; }

  ScalarType secondDerivativeImpl(const JType& /* J */) const { return 1; };
};

template <typename ST>
struct VF2T
{
  using ScalarType = ST;
  using JType      = ScalarType;

  ScalarType storedEnergyImpl(const JType& J) const { return 0.25 * (std::pow(J - 1, 2) + std::pow(std::log(J), 2)); };

  ScalarType firstDerivativeImpl(const JType& J) const { return 0.5 * (J - 1 + 1 / J * log(J)); }

  ScalarType secondDerivativeImpl(const JType& J) const {
    auto Jsq = std::pow(J, 2);
    return 1 / (2 * Jsq) * (1 + Jsq - log(J));
  }
};

template <typename ST>
struct VF3T
{
  using ScalarType = ST;
  using JType      = ScalarType;

  ScalarType storedEnergyImpl(const JType& J) const { return 0.5 * std::pow(log(J), 2); };

  ScalarType firstDerivativeImpl(const JType& J) const { return 1 / J * log(J); }

  ScalarType secondDerivativeImpl(const JType& J) const {
    auto Jsq = std::pow(J, 2);
    return 1 / Jsq * (1 - log(J));
  }
};

template <typename ST>
struct VF4T
{
  using ScalarType = ST;
  using JType      = ScalarType;

  ScalarType storedEnergyImpl(const JType& J) const {
    return (1 / std::pow(beta_, 2)) * ((1 / std::pow(J, beta_)) - 1 + beta_ * std::log(J));
  };

  ScalarType firstDerivativeImpl(const JType& J) const { return (1 / beta_) * ((1 / J) - (1 / (pow(J, 1 + beta_)))); }

  ScalarType secondDerivativeImpl(const JType& J) const {
    return (1 / beta_) * ((1 / std::pow(J, 2 + beta_)) * (1 + beta_ - std::pow(J, beta_)));
  }

  explicit VF4T(double beta)
      : beta_(beta) {}

private:
  double beta_;
};

template <typename ST>
struct VF5T
{
  using ScalarType = ST;
  using JType      = ScalarType;

  ScalarType storedEnergyImpl(const JType& J) const { return 0.25 * (std::pow(J, 2) - 1 - 2 * std::log(J)); };

  ScalarType firstDerivativeImpl(const JType& J) const { return 0.5 * (J - (1 / J)); }

  ScalarType secondDerivativeImpl(const JType& J) const { return 0.5 * (1 + (1 / std::pow(J, 2))); }
};

template <typename ST>
struct VF6T
{
  using ScalarType = ST;
  using JType      = ScalarType;

  ScalarType storedEnergyImpl(const JType& J) const { return J - std::log(J) - 1; };

  ScalarType firstDerivativeImpl(const JType& J) const { return 1 - (1 / J); }

  ScalarType secondDerivativeImpl(const JType& J) const { return 1 / std::pow(J, 2); }
};

using VF1 = VF1T<double>;
using VF2 = VF2T<double>;
using VF3 = VF3T<double>;
using VF4 = VF4T<double>;
using VF5 = VF5T<double>;
using VF6 = VF6T<double>;


// Default implementation (no volumetric part)
template <typename ST>
struct VF0T
{
  using ScalarType = ST;
  using JType      = ScalarType;

  ScalarType storedEnergyImpl(const JType& /* J */) const { return 0; };

  ScalarType firstDerivativeImpl(const JType& /* J */) const { return 0; }

  ScalarType secondDerivativeImpl(const JType& /* J */) const { return 0; };
};

template <typename ST>
struct VolumetricPart<VF0T<ST>>
{
  using ScalarType = ST;
  using JType      = typename VF0T<ST>::JType;

  ScalarType storedEnergy(const JType& /* J */) const { return 0; };

  ScalarType firstDerivative(const JType& /* J */) const { return 0; }

  ScalarType secondDerivative(const JType& /* J */) const { return 0; };
};
using NoVolumetricPart = VolumetricPart<VF0T<double>>;

} // namespace Ikarus
