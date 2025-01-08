// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file factory.hh
 * \brief Header file for containing helper functions to construct various hyperelastic material models.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/ogden.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/neohooke.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/volumetric/volumetricfunctions.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>

namespace Ikarus::Materials {

/**
 * \brief A helper function to create a hyperelastic material model, where the deviatoric part is according to the Ogden
 * model.
 *
 * \tparam n Number of Ogden material parameters.
 * \tparam PrincipalStretchTag Type of the principal stretches, either total or deviatoric.
 * \tparam VolumetricFunction Type of the volumetric function.
 *
 * \param mu The shear parameters (mu_i).
 * \param alpha The (exponential) parameters (alpha_i).
 * \param K Bulk modulus (or) Lamé's first parameter.
 * \param vf The volumetric function.
 *
 * \return A hyperelastic material model.
 */
template <int n, PrincipalStretchTag tag, typename VolumetricFunction = VF0T<double>>
inline auto makeOgden(const typename Ogden<n, tag>::MaterialParameters& mu,
                      const typename Ogden<n, tag>::MaterialExponents alpha, double K = 0.0,
                      const VolumetricFunction& vf = VolumetricFunction{}) {
  auto ogPre = Ogden<n, tag>(mu, alpha);
  auto dev   = Deviatoric(ogPre);
  auto vol   = Volumetric(K, vf);

  return Hyperelastic(dev, vol);
}

} // namespace Ikarus::Materials
