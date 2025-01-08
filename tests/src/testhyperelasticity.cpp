// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhyperelasticity.hh"

#include <dune/common/test/testsuite.hh>

#include <ikarus/utils/init.hh>

using Dune::TestSuite;

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testVolumetricFunctions());
  t.subTest(recoverNeoHookeTest());
  t.subTest(testMaterialResults());

  return t.exit();
}
