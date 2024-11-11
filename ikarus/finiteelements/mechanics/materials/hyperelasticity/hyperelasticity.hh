// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file materials.hh
 * \brief Header file for material models in Ikarus finite element mechanics.
 * \ingroup  mechanics
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/deviatoric.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/hyperelastic.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/modified/blatzko.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/modified/ogden.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/regularized/ogden.hh>


namespace Ikarus {
auto makeBlatzKo(ShearModulus mu) {
  auto bk  = BlatzKo(mu);
  auto dev = Deviatoric<decltype(bk)>(bk);

  return Hyperelastic(dev);
}

MAKE_ENUM(RegularizedTag, regularized, modified);

template <int n, RegularizedTag tag, typename VolumetricFunction = VF0T<double>>
auto makeOgden(const typename RegOgden<n>::MaterialParameters& mu, const typename RegOgden<n>::OgdenParameters og,
               BulkModulus K = {0.0}, const VolumetricFunction& vf = VolumetricFunction{}) {
  if constexpr (tag == RegularizedTag::regularized)
    return makeRegularizedOgden<n>(mu, og, K, vf);
  else
    return makeModifiedOgden<n>(mu, og, K, vf);
}

/**
 * \brief This returns the compressible version of the ogden material using the regularized stretches
 */
template <int n, typename VolumetricFunction = VF0T<double>>
auto makeRegularizedOgden(const typename RegOgden<n>::MaterialParameters& mu,
                          const typename RegOgden<n>::OgdenParameters og, BulkModulus K = {0.0},
                          const VolumetricFunction& vf = VolumetricFunction{}) {
  auto ogPre = RegOgden<n>(mu, og);
  auto dev   = Deviatoric<decltype(ogPre)>(ogPre);
  auto vol   = Volumetric(K, vf);

  return Hyperelastic(dev, vol);
}

/**
 * \brief This returns the compressible version of the ogden material using the unmodified stretches and a addtional
 * component in the energy of  -mu * ln(J) to ensure a stress free reference configuration
 */
template <int n, typename VolumetricFunction = VF0T<double>>
auto makeModifiedOgden(const typename ModOgden<n>::MaterialParameters& mu,
                       const typename ModOgden<n>::OgdenParameters og, BulkModulus K = {0.0},
                       const VolumetricFunction& vf = VolumetricFunction{}) {
  auto ogPre = ModOgden<n>(mu, og);
  auto dev   = Deviatoric<decltype(ogPre)>(ogPre);
  auto vol   = Volumetric(K, vf);

  return Hyperelastic(dev, vol);
}

} // namespace Ikarus
