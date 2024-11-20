// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
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
  auto bk = makeBlatzKo(ShearModulus{mu});

  std::array<double, 3> mu_og    = {2.0 * mu / 3.0, mu / 6.0, mu / 6.0};
  std::array<double, 3> alpha_og = {1.23, 0.59, 0.18};
  auto ogdenTotal                = makeOgden<3, PrincipalStretchTag::total>(mu_og, alpha_og, {Lambda}, VF3{});
  auto ogdenDevi                 = makeOgden<3, PrincipalStretchTag::deviatoric>(mu_og, alpha_og, {K}, VF3{});
  auto mr                        = makeMooneyRivlin({mu / 2.0, mu / 2.0}, {K}, VF3{});
  auto yeoh                      = makeYeoh({2.0 * mu / 3.0, mu / 6.0, mu / 6.0}, {K}, VF3{});
  auto polynomial = makeInvariantBased<5>({mu / 6.0, mu / 8.0, mu / 8.0, mu / 12.0, mu / 2.0}, {0, 1, 3, 6, 8},
                                          {2, 4, 10, 1, 7}, {K}, VF3{});
  auto ab         = makeArrudaBoyce({mu, 0.85}, {K}, VF3{});
  auto gent       = makeGent({mu, 2.5}, {K}, VF3{});

  t.subTest(checkMaterialByAutoDiff(nh));
  t.subTest(checkMaterialByAutoDiff(bk));
  t.subTest(checkMaterialByAutoDiff(ogdenTotal));
  t.subTest(checkMaterialByAutoDiff(ogdenDevi));
  t.subTest(checkMaterialByAutoDiff(mr));
  t.subTest(checkMaterialByAutoDiff(yeoh));
  t.subTest(checkMaterialByAutoDiff(polynomial));
  t.subTest(checkMaterialByAutoDiff(ab));
  t.subTest(checkMaterialByAutoDiff(gent));

  return t.exit();
}
