// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file materials.hh
 * \brief Header file for material models in Ikarus finite element mechanics.
 * \ingroup  mechanics
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/blatzko.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/deviatoric.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/hyperelastic.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/modogden.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/ogden.hh>

namespace Ikarus {
auto makeBlatzKo(ShearModulus mu) {
  auto bk  = BlatzKo(mu);
  auto dev = DeviatoricPart<decltype(bk)>(bk);

  return Hyperelastic(dev);
}

/**
 * \brief This returns the compressible version of the ogden material using the modified stretches
 */
template <int n, typename VolumetricFunction = VF0T<double>>
auto makeOgden(const typename Ogden<n>::MaterialParameters& mu, const typename Ogden<n>::OgdenParameters og,
               BulkModulus K = {0.0}, const VolumetricFunction& vf = VolumetricFunction{}) {
  auto ogPre = Ogden<n>(mu, og);
  auto dev   = DeviatoricPart<decltype(ogPre)>(ogPre);
  auto vol   = VolumetricPart(K, vf);

  return Hyperelastic(dev, vol);
}

/**
 * \brief This returns the compressible version of the ogden material using the unmodified stretches and a addtional
 * component in the energy of  -mu * ln(J) to ensure a stress free reference configuration
 */
template <int n, typename VolumetricFunction = VF0T<double>>
auto makeModifiedOgden(const typename ModifiedOgden<n>::MaterialParameters& mu,
                       const typename ModifiedOgden<n>::OgdenParameters og, BulkModulus K = {0.0},
                       const VolumetricFunction& vf = VolumetricFunction{}) {
  auto ogPre = ModifiedOgden<n>(mu, og);
  auto dev   = DeviatoricPart<decltype(ogPre)>(ogPre);
  auto vol   = VolumetricPart(K, vf);

  return Hyperelastic(dev, vol);
}

} // namespace Ikarus
