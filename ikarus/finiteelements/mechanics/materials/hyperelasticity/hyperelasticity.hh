// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file hyperelasticity.hh
 * \brief Header file for material models in Ikarus finite element mechanics.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/deviatoric.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/hyperelastic.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/modified/blatzko.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/modified/ogden.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/regularized/ogden.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>

namespace Ikarus::Materials {
auto makeBlatzKo(ShearModulus mu) {
  auto bk  = BlatzKo(mu);
  auto dev = Deviatoric<decltype(bk)>(bk);

  return Hyperelastic(dev);
}


template <int n, StretchTag tag, typename VolumetricFunction = VF0T<double>>
auto makeOgden(const typename RegOgden<n>::MaterialParameters& mu, const typename RegOgden<n>::OgdenParameters og,
               BulkModulus K = {0.0}, const VolumetricFunction& vf = VolumetricFunction{}) {
  using OgdenType = std::conditional_t<tag == StretchTag::deviatoric, RegOgden<n>, ModOgden<n>>;
  auto ogPre      = OgdenType(mu, og);
  auto dev        = Deviatoric<decltype(ogPre)>(ogPre);
  auto vol        = Volumetric(K, vf);

  return Hyperelastic(dev, vol);
}

} // namespace Ikarus::Materials
