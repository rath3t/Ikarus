// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <dune/common/test/testsuite.hh>
using Dune::TestSuite;
#include <dune/common/float_cmp.hh>

#include <ikarus/utils/init.hh>
#include <ikarus/utils/tensorproductquadrule.hh>

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  for (int i = 0; i < 5; ++i) {
    const auto& twoDRule = Dune::QuadratureRules<double, 2>::rule(Dune::GeometryTypes::quadrilateral, i);
    const auto& oneDRule = Dune::QuadratureRules<double, 1>::rule(Dune::GeometryTypes::line, 1);

    const auto rule = Ikarus::tensorProductQuadrature(twoDRule, oneDRule);
    double vol      = 0;
    for (int gpIndex = 0; auto& gp : rule) {
      vol += gp.weight();
      gpIndex++;
    }
    t.check(Dune::FloatCmp::eq(vol, 1.0));
  }
}
