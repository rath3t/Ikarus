// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file volumetricfunctions.hh
 * \brief Implementation of the volumetric part of a hyperelastic material.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/concepts.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/volumetric/interface.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

// Default implementation (no volumetric part)
template <typename ST>
struct VF0T
{
  using ScalarType = ST;
  using JType      = ScalarType;

  ScalarType storedEnergyImpl(const JType& /* J */) const { return 0; };

  ScalarType firstDerivativeImpl(const JType& /* J */) const { return 0; }

  ScalarType secondDerivativeImpl(const JType& /* J */) const { return 0; };

  template <typename STO>
  auto rebind() const {
    return VF0T<STO>();
  }

  [[nodiscard]] constexpr static std::string name() noexcept { return "None"; }
};
/**
 * \brief Volumetric function No. 1 found in \cite hartmann_polyconvexity_2003 Tab. 4
 * \details $U(J) = \frac{1}{2}(J - 1)^2$
 * \tparam ST ScalarType
 * \ingroup materials
 */
template <typename ST>
struct VF1T
{
  using ScalarType = ST;
  using JType      = ScalarType;

  ScalarType storedEnergyImpl(const JType& J) const { return 0.5 * pow(J - 1, 2); };

  ScalarType firstDerivativeImpl(const JType& J) const { return J - 1; }

  ScalarType secondDerivativeImpl(const JType& /* J */) const { return 1; };

  template <typename STO>
  auto rebind() const {
    return VF1T<STO>();
  }

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 1"; }
};

/**
 * \brief Volumetric function No. 2 found in \cite hartmann_polyconvexity_2003 Tab. 4
 * \details $U(J) = \frac{1}{4}\left((J - 1)^2 + (\ln J )^2 \right)$ according to Simo and Taylor 1982
 * \tparam ST ScalarType
 * \ingroup materials
 */
template <typename ST>
struct VF2T
{
  using ScalarType = ST;
  using JType      = ScalarType;

  ScalarType storedEnergyImpl(const JType& J) const { return 0.25 * (pow(J - 1, 2) + pow(log(J), 2)); };

  ScalarType firstDerivativeImpl(const JType& J) const { return 0.5 * (J - 1 + 1 / J * log(J)); }

  ScalarType secondDerivativeImpl(const JType& J) const {
    auto Jsq = pow(J, 2);
    return 1 / (2 * Jsq) * (1 + Jsq - log(J));
  }

  template <typename STO>
  auto rebind() const {
    return VF2T<STO>();
  }

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 2"; }
};

template <typename ST>
struct VF3T
{
  using ScalarType = ST;
  using JType      = ScalarType;

  ScalarType storedEnergyImpl(const JType& J) const { return 0.5 * pow(log(J), 2); };

  ScalarType firstDerivativeImpl(const JType& J) const { return 1 / J * log(J); }

  ScalarType secondDerivativeImpl(const JType& J) const { return 1 / pow(J, 2) * (1 - log(J)); }

  template <typename STO>
  auto rebind() const {
    return VF3T<STO>();
  }

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 3"; }
};

template <typename ST>
struct VF4T
{
  using ScalarType = ST;
  using JType      = ScalarType;

  ScalarType storedEnergyImpl(const JType& J) const {
    return (1 / pow(beta_, 2)) * ((1 / pow(J, beta_)) - 1 + beta_ * log(J));
  };

  ScalarType firstDerivativeImpl(const JType& J) const { return (1 / beta_) * ((1 / J) - (1 / (pow(J, 1 + beta_)))); }

  ScalarType secondDerivativeImpl(const JType& J) const {
    return (1 / beta_) * ((1 / pow(J, 2 + beta_)) * (1 + beta_ - pow(J, beta_)));
  }

  explicit VF4T(double beta)
      : beta_(beta) {}

  template <typename STO>
  auto rebind() const {
    return VF4T<STO>(beta_);
  }

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 4"; }

  double beta() const { return beta_; }

private:
  double beta_;
};

template <typename ST>
struct VF5T
{
  using ScalarType = ST;
  using JType      = ScalarType;

  ScalarType storedEnergyImpl(const JType& J) const { return 0.25 * (pow(J, 2) - 1 - 2 * log(J)); };

  ScalarType firstDerivativeImpl(const JType& J) const { return 0.5 * (J - (1 / J)); }

  ScalarType secondDerivativeImpl(const JType& J) const { return 0.5 * (1 + (1 / pow(J, 2))); }

  template <typename STO>
  auto rebind() const {
    return VF5T<STO>();
  }

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 5"; }
};

template <typename ST>
struct VF6T
{
  using ScalarType = ST;
  using JType      = ScalarType;

  ScalarType storedEnergyImpl(const JType& J) const { return J - log(J) - 1; };

  ScalarType firstDerivativeImpl(const JType& J) const { return 1 - (1 / J); }

  ScalarType secondDerivativeImpl(const JType& J) const { return 1 / pow(J, 2); }

  template <typename STO>
  auto rebind() const {
    return VF6T<STO>();
  }

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 6"; }
};

using NoVolumetricPart = Volumetric<VF0T<double>>;
using VF1              = VF1T<double>;
using VF2              = VF2T<double>;
using VF3              = VF3T<double>;
using VF4              = VF4T<double>;
using VF5              = VF5T<double>;
using VF6              = VF6T<double>;

} // namespace Ikarus::Materials
