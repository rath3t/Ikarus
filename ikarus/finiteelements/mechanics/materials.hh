// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file materials.hh
 * \brief Header file for material models in Ikarus finite element mechanics.
 * \ingroup  mechanics
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/blatzko.hh>
#include <ikarus/finiteelements/mechanics/materials/compressibleogden.hh>
#include <ikarus/finiteelements/mechanics/materials/deviatoric.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic.hh>
#include <ikarus/finiteelements/mechanics/materials/incompressibleogden.hh>
#include <ikarus/finiteelements/mechanics/materials/linearelasticity.hh>
#include <ikarus/finiteelements/mechanics/materials/neohooke.hh>
#include <ikarus/finiteelements/mechanics/materials/svk.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/finiteelements/mechanics/materials/vanishingstrain.hh>
#include <ikarus/finiteelements/mechanics/materials/vanishingstress.hh>

namespace Ikarus {
auto makeBlatzKo(ShearModulus mu) {
  auto bk  = BlatzKo(mu);
  auto dev = DeviatoricPart<decltype(bk), false>(bk);

  return Hyperelastic(dev);
}

template <int n, typename VolumetricFunction = VF0T<double>>
auto makeIncompressibleOgden(const typename IncompressibleOgden<n>::MaterialParameters& mu,
                             const typename IncompressibleOgden<n>::OgdenParameters og, BulkModulus K = {0.0},
                             const VolumetricFunction& vf = VolumetricFunction{}) {
  auto ogPre = IncompressibleOgden<n>(mu, og);
  auto dev = DeviatoricPart<decltype(ogPre), true>(ogPre);
  auto vol = VolumetricPart(K, vf);

  return Hyperelastic(dev, vol);

}

template <int n, typename VolumetricFunction = VF0T<double>>
auto makeCompressibleOgden(const typename IncompressibleOgden<n>::MaterialParameters& mu,
                             const typename IncompressibleOgden<n>::OgdenParameters og, BulkModulus K = {0.0},
                             const VolumetricFunction& vf = VolumetricFunction{}) {
  auto ogPre = CompressibleOgden<n>(mu, og);
  auto dev = DeviatoricPart<decltype(ogPre), false>(ogPre);
  auto vol = VolumetricPart(K, vf);

  return Hyperelastic(dev, vol);

}

} // namespace Ikarus
