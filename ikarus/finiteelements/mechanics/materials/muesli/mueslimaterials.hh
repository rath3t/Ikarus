// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file materials.hh
 * \brief Header file for material models in Ikarus finite element mechanics.
 * \ingroup  mechanics
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/muesli/muesliefinite.hh>
#include <ikarus/finiteelements/mechanics/materials/muesli/muesliesmall.hh>
#include <ikarus/finiteelements/mechanics/materials/muesli/mueslihelpers.hh>

namespace Ikarus::Materials::Muesli {

// Alias for Muesli materials

using LinearElasticity            = muesli::elasticIsotropicMaterial;
using LinearAnisotropicElasticity = muesli::elasticAnisotropicMaterial;
using LinearOrthotropicElasticity = muesli::elasticOrthotropicMaterial;

using StVenantKirchhoff = muesli::svkMaterial;
using NeoHooke          = muesli::neohookeanMaterial;
using MooneyRivlin      = muesli::mooneyMaterial;
using Yeoh              = muesli::yeohMaterial;
using ArrudaBoyce       = muesli::arrudaboyceMaterial;

template <typename MPT>
auto makeLinearElasticity(const MPT& mpt) {
  auto muesliParameters = propertiesFromIkarusMaterialParameters(mpt);
  return SmallStrain<>(muesliParameters);
}

template <typename MPT>
auto makeNeoHooke(const MPT& mpt, bool useDeviatoricStretches = false) {
  auto muesliParameters = propertiesFromIkarusMaterialParameters(mpt);
  if (useDeviatoricStretches)
    addRegularizedTag(muesliParameters);
  return FiniteStrain<Muesli::NeoHooke>(muesliParameters);
}

template <typename MPT>
auto makeSVK(const MPT& mpt) {
  auto muesliParameters = propertiesFromIkarusMaterialParameters(mpt);
  return FiniteStrain<Muesli::StVenantKirchhoff>(muesliParameters);
}

inline auto makeArrudaBoyce(double C1, double lambda_m, double K, bool compressible = true) {
  auto muesliParameters = muesli::materialProperties{};
  muesliParameters.insert({"c1", C1});
  muesliParameters.insert({"lambdam", lambda_m});
  muesliParameters.insert({"bulk", K});
  if (compressible)
    addCompressibleTag(muesliParameters);
  return FiniteStrain<Muesli::ArrudaBoyce>(muesliParameters);
}

inline auto makeYeoh(std::array<double, 3> C, double K, bool compressible = true) {
  auto muesliParameters = muesli::materialProperties{};
  muesliParameters.insert({"c1", C[0]});
  muesliParameters.insert({"c2", C[1]});
  muesliParameters.insert({"c3", C[2]});
  muesliParameters.insert({"bulk", K});
  if (compressible)
    addCompressibleTag(muesliParameters);
  return FiniteStrain<Muesli::Yeoh>(muesliParameters);
}

inline auto makeMooneyRivlin(std::array<double, 3> alpha, bool incompressible = false) {
  auto muesliParameters = muesli::materialProperties{};
  muesliParameters.insert({"alpha0", alpha[0]});
  muesliParameters.insert({"alpha1", alpha[1]});
  muesliParameters.insert({"alpha2", alpha[2]});
  if (incompressible)
    addIncompressibleTag(muesliParameters);
  return FiniteStrain<Muesli::MooneyRivlin>(muesliParameters);
}

} // namespace Ikarus::Materials::Muesli