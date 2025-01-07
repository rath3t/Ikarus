// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <dune/common/test/testsuite.hh>

#include <ikarus/utils/init.hh>
#include <ikarus/utils/pythonautodiffdefinitions.hh>

using Dune::TestSuite;

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  autodiff::real valDual = 7.0;
  valDual[1]             = 1;
  Python::start();
  auto pyLambda = Python::Conversion<autodiff::Real<1, double>>::toPy(valDual);

  autodiff::real valExpected;
  Python::Conversion<autodiff::Real<1, double>>::toC(pyLambda, valExpected);

  t.check(valDual == valExpected);

  return t.exit();
}
