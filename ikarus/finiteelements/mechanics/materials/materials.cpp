// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file materials.cpp
 * \brief Explicit instantiations for material templates.
 */

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/factory.hh>
#include <ikarus/finiteelements/mechanics/materials/linearelasticity.hh>
#include <ikarus/finiteelements/mechanics/materials/svk.hh>

namespace Ikarus::Materials {

template struct LinearElasticityT<double>;
template struct NeoHookeT<double>;
template struct StVenantKirchhoffT<double>;
template struct Hyperelastic<Deviatoric<BlatzKoT<double>>>;

// Here only the deviatoric parts can be precompiled
template struct GentT<double>;
template struct InvariantBasedT<double, 2>; // Money-Rivlin
template struct InvariantBasedT<double, 3>; // Yeoh
template struct ArrudaBoyceT<double>;

} // namespace Ikarus::Materials
