// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file factory.hh
 * \brief Header file for hyperelastic material models in Ikarus finite element mechanics.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/blatzko.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/ogden.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/neohooke.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/volumetric/volumetricfunctions.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>

namespace Ikarus::Materials {

inline auto makeBlatzKo(ShearModulus mu) {
  auto bk  = BlatzKo(mu);
  auto dev = Deviatoric<decltype(bk)>(bk);

  return Hyperelastic(dev);
}

template <int n, PrincipalStretchTag tag, typename VolumetricFunction = VF0T<double>>
inline auto makeOgden(const typename Ogden<n, tag>::MaterialParameters& mu,
                      const typename Ogden<n, tag>::OgdenParameters og, BulkModulus K = {0.0},
                      const VolumetricFunction& vf = VolumetricFunction{}) {
  auto ogPre = Ogden<n, tag>(mu, og);
  auto dev   = Deviatoric(ogPre);
  auto vol   = Volumetric(K, vf);

  return Hyperelastic(dev, vol);
}

} // namespace Ikarus::Materials
