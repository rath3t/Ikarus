// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file factory.hh
 * \brief Header file for containing helper functions to construct various hyperelastic material models.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/arrudaboyce.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/blatzko.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/gent.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/invariantbased.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/ogden.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/neohooke.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/volumetric/volumetricfunctions.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>

namespace Ikarus::Materials {

/**
 * \brief A helper function to create a hyperelastic material model, where the deviatoric part is according to the
 * Blatz-Ko model.
 *
 *\param mu The shear modulus.
 *
 *\return A hyperelastic material model.
 */
inline auto makeBlatzKo(double mu) {
  auto bk  = BlatzKo(mu);
  auto dev = Deviatoric<decltype(bk)>(bk);

  return Hyperelastic(dev);
}

/**
 * \brief A helper function to create a hyperelastic material model, where the deviatoric part is according to the Ogden
 * model.
 *
 * \tparam n Number of Ogden material parameters.
 * \tparam PrincipalStretchTag Type of the principal stretches, either total or deviatoric.
 * \tparam VolumetricFunction Type of the volumetric function.
 *
 * \param mu The shear parameters (mu_i).
 * \param og The (exponential) parameters (alpha_i).
 * \param K The bulk modulus.
 * \param vf The volumetric function.
 *
 * \return A hyperelastic material model.
 */
template <int n, PrincipalStretchTag tag, typename VolumetricFunction = VF0T<double>>
inline auto makeOgden(const typename Ogden<n, tag>::MaterialParameters& mu,
                      const typename Ogden<n, tag>::OgdenParameters og, double K = 0.0,
                      const VolumetricFunction& vf = VolumetricFunction{}) {
  auto ogPre = Ogden<n, tag>(mu, og);
  auto dev   = Deviatoric(ogPre);
  auto vol   = Volumetric(K, vf);

  return Hyperelastic(dev, vol);
}

/**
 * \brief A helper function to create a hyperelastic material model, where the deviatoric part is according to the
 * general invariant-based model.
 *
 * \tparam n Number of material parameters.
 * \tparam VolumetricFunction Type of the volumetric function.
 *
 * \param mu The shear parameters (mu_i).
 * \param pex The exponents related to the first invariant.
 * \param qex The exponents related to the second invariant.
 * \param K The bulk modulus.
 * \param vf The volumetric function.
 *
 * \return A hyperelastic material model.
 */
template <int n, typename VolumetricFunction = VF0T<double>>
inline auto makeInvariantBased(const typename InvariantBased<n>::MaterialParameters& mu,
                               const typename InvariantBased<n>::Exponents pex,
                               const typename InvariantBased<n>::Exponents qex, double K = 0.0,
                               const VolumetricFunction& vf = VolumetricFunction{}) {
  auto invariantBasedPre = InvariantBased<n>(pex, qex, mu);
  auto dev               = Deviatoric(invariantBasedPre);
  auto vol               = Volumetric(K, vf);

  return Hyperelastic(dev, vol);
}

/**
 * \brief A helper function to create a hyperelastic material model, where the deviatoric part is according to the
 * Mooney-Rivlin model (InvariantBased model with n = 2).
 *
 * \tparam VolumetricFunction Type of the volumetric function.
 *
 * \param mu The shear parameters (mu_i).
 * \param K The bulk modulus.
 * \param vf The volumetric function.
 *
 * \return A hyperelastic material model.
 */
template <typename VolumetricFunction = VF0T<double>>
inline auto makeMooneyRivlin(const typename InvariantBased<2>::MaterialParameters& mu, double K = 0.0,
                             const VolumetricFunction& vf = VolumetricFunction{}) {
  typename InvariantBased<2>::Exponents pex = {1, 0};
  typename InvariantBased<2>::Exponents qex = {0, 1};

  return makeInvariantBased<2, VolumetricFunction>(mu, pex, qex, K, vf);
}

/**
 * \brief A helper function to create a hyperelastic material model, where the deviatoric part is according to the
 * Yeoh model (InvariantBased model with n = 3).
 *
 * \tparam VolumetricFunction Type of the volumetric function.
 *
 * \param mu The shear parameters (mu_i).
 * \param K The bulk modulus.
 * \param vf The volumetric function.
 *
 * \return A hyperelastic material model.
 */
template <typename VolumetricFunction = VF0T<double>>
inline auto makeYeoh(const typename InvariantBased<3>::MaterialParameters& mu, double K = 0.0,
                     const VolumetricFunction& vf = VolumetricFunction{}) {
  typename InvariantBased<3>::Exponents pex = {1, 2, 3};
  typename InvariantBased<3>::Exponents qex = {0, 0, 0};

  return makeInvariantBased<3, VolumetricFunction>(mu, pex, qex, K, vf);
}

/**
 * \brief A helper function to create a hyperelastic material model, where the deviatoric part is according to the
 * Arruda-Boyce model.
 *
 * \tparam VolumetricFunction Type of the volumetric function.
 *
 * \param matPar The Arruda-Boyce material parameters (C and lambdaM).
 * \param K The bulk modulus.
 * \param vf The volumetric function.
 *
 * \return A hyperelastic material model.
 */
template <typename VolumetricFunction = VF0T<double>>
inline auto makeArrudaBoyce(const ArrudaBoyceMatParameters& matPar, double K = 0.0,
                            const VolumetricFunction& vf = VolumetricFunction{}) {
  auto abPre = ArrudaBoyce(matPar);
  auto dev   = Deviatoric(abPre);
  auto vol   = Volumetric(K, vf);

  return Hyperelastic(dev, vol);
}

/**
 * \brief A helper function to create a hyperelastic material model, where the deviatoric part is according to the
 * Gent model.
 *
 * \tparam VolumetricFunction Type of the volumetric function.
 *
 * \param matPar The Gent material parameters (mu and Jm).
 * \param K The bulk modulus.
 * \param vf The volumetric function.
 *
 * \return A hyperelastic material model.
 */
template <typename VolumetricFunction = VF0T<double>>
inline auto makeGent(const GentMatParameters& matPar, double K = 0.0,
                     const VolumetricFunction& vf = VolumetricFunction{}) {
  auto gentPre = Gent(matPar);
  auto dev     = Deviatoric(gentPre);
  auto vol     = Volumetric(K, vf);

  return Hyperelastic(dev, vol);
}

} // namespace Ikarus::Materials
