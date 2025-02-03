// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file materials.hh
 * \brief Header file for material models in Ikarus finite element mechanics.
 * \ingroup  mechanics
 */

#pragma once

#if ENABLE_MUESLI

  #include <ikarus/finiteelements/mechanics/materials/muesli/mueslifinite.hh>
  #include <ikarus/finiteelements/mechanics/materials/muesli/mueslihelpers.hh>
  #include <ikarus/finiteelements/mechanics/materials/muesli/mueslismall.hh>

namespace Ikarus::Materials {

// Alias for small strain Muesli materials

/**
 * \brief Constructs a muesli linear isotropic elastic material model based on muesli::elasticIsotropicMaterial.
 *
 * \tparam MPT the type of the Ikarus material parameters.
 * \param mpt the Ikarus material parameters.
 * \return auto the constructed material properties.
 */
template <Concepts::MPTuple MPT>
inline auto makeMuesliLinearElasticity(const MPT& mpt) {
  auto muesliParameters = propertiesFromIkarusMaterialParameters(mpt);
  return SmallStrain<muesli::elasticIsotropicMaterial>(muesliParameters);
}

/**
 * \brief Constructs a muesli NeoHooke material model based on muesli::neohookeanMaterial.
 *
 * \tparam MPT the type of the Ikarus material parameters.
 * \param mpt the Ikarus material parameters.
 * \param useDeviatoricStretches tells the material to use deviatoric principal stretches.
 * \return auto the constructed material.
 */
template <Concepts::MPTuple MPT>
inline auto makeMuesliNeoHooke(const MPT& mpt, bool useDeviatoricStretches = false) {
  auto muesliParameters = propertiesFromIkarusMaterialParameters(mpt);
  if (useDeviatoricStretches)
    addTag(muesliParameters, "subtype regularized");
  return FiniteStrain<muesli::neohookeanMaterial>(muesliParameters);
}

/**
 * \brief Constructs a muesli StVenantKirchhoff material model based on muesli::svkMaterial.
 *
 * \tparam MPT the type of the Ikarus material parameters.
 * \param mpt the Ikarus material parameters.
 * \return auto the constructed material.
 */
template <Concepts::MPTuple MPT>
inline auto makeMuesliSVK(const MPT& mpt) {
  auto muesliParameters = propertiesFromIkarusMaterialParameters(mpt);
  return FiniteStrain<muesli::svkMaterial>(muesliParameters);
}

/**
 * \brief Constructs a muesli ArrudaBoyce material model based on muesli::arrudaboyceMaterial.
 *
 * \param C1 C1 parameter.
 * \param lambda_m lambda_m parameter.
 * \param K the bulk modulus.
 * \param compressible tells the material to use the compressible version of ArrudaBoyce.
 * \return auto the constructed material.
 */
inline auto makeMuesliArrudaBoyce(double C1, double lambda_m, double K, bool compressible = true) {
  auto muesliParameters = muesli::materialProperties{};
  muesliParameters.insert({"c1", C1});
  muesliParameters.insert({"lambdam", lambda_m});
  muesliParameters.insert({"bulk", K});
  if (compressible)
    addTag(muesliParameters, "compressible");
  return FiniteStrain<muesli::arrudaboyceMaterial>(muesliParameters);
}

/**
 * \brief Constructs a muesli Yeoh material model based on muesli::yeohMaterial.
 *
 * \param C c1, c2, c3 parameters in a std::array.
 * \param K the bzlk modolus.
 * \param compressible tells the material to use the compressible version of Yeoh.
 * \return auto the constructed material.
 */
inline auto makeMuesliYeoh(std::array<double, 3> C, double K, bool compressible = true) {
  auto muesliParameters = muesli::materialProperties{};
  muesliParameters.insert({"c1", C[0]});
  muesliParameters.insert({"c2", C[1]});
  muesliParameters.insert({"c3", C[2]});
  muesliParameters.insert({"bulk", K});
  if (compressible)
    addTag(muesliParameters, "compressible");
  return FiniteStrain<muesli::yeohMaterial>(muesliParameters);
}

/**
 * \brief Constructs a muesli MooneyRivlin material model based on muesli::mooneyMaterial.
 *
 * \param alpha alpha0, alpha1, alpha2 parameters in a std::array, where alpha0 is the bulk modulus.
 * \param incompressible tells the material to use the incompressible version of MooneyRivlin.
 * \return auto the constructed material.
 */
inline auto makeMooneyRivlin(std::array<double, 3> alpha, bool incompressible = false) {
  auto muesliParameters = muesli::materialProperties{};
  muesliParameters.insert({"alpha0", alpha[0]});
  muesliParameters.insert({"alpha1", alpha[1]});
  muesliParameters.insert({"alpha2", alpha[2]});
  if (incompressible)
    addTag(muesliParameters, "incompressible");
  return FiniteStrain<muesli::mooneyMaterial>(muesliParameters);
}

} // namespace Ikarus::Materials

#else
  #error Muesli materials depends on the Muesli library, which is not included
#endif