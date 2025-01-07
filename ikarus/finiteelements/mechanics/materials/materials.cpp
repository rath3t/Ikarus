// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file materials.cpp
 * \brief Explicit instantiations for material templates.
 */

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/neohooke.hh>
#include <ikarus/finiteelements/mechanics/materials/linearelasticity.hh>
#include <ikarus/finiteelements/mechanics/materials/svk.hh>

namespace Ikarus::Materials {

template struct LinearElasticityT<double>;
template struct NeoHookeT<double>;
template struct StVenantKirchhoffT<double>;

} // namespace Ikarus::Materials
