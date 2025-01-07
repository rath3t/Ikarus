// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "checkmaterialbyautodiff.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Eigenvalues>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/init.hh>

using namespace Ikarus;
using namespace Ikarus::Materials;
using Dune::TestSuite;

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  // instantiate material models
  double Emod     = 1000;
  double nu       = 0.25; // Blatz Ko assumes nu = 0.25
  auto matPar     = YoungsModulusAndPoissonsRatio{.emodul = Emod, .nu = nu};
  auto mu         = convertLameConstants(matPar).toShearModulus();
  auto K          = convertLameConstants(matPar).toBulkModulus();
  auto Lambda     = convertLameConstants(matPar).toLamesFirstParameter();
  auto lameMatPar = toLamesFirstParameterAndShearModulus(matPar);

  auto nh = NeoHooke(lameMatPar);

  std::array<double, 3> mu_og    = {2.0 * mu / 3.0, mu / 6.0, mu / 6.0};
  std::array<double, 3> alpha_og = {1.23, 0.59, 0.18};
  auto ogdenTotal                = makeOgden<3, PrincipalStretchTag::total>(mu_og, alpha_og, Lambda, VF3{});
  auto ogdenDevi                 = makeOgden<3, PrincipalStretchTag::deviatoric>(mu_og, alpha_og, K, VF3{});

  t.subTest(checkMaterialByAutoDiff(nh));
  t.subTest(checkMaterialByAutoDiff(ogdenTotal));
  t.subTest(checkMaterialByAutoDiff(ogdenDevi));

  return t.exit();
}
