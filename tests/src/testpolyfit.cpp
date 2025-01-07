// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <matplot/matplot.h>

#include <dune/common/float_cmp.hh>
#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/utils/init.hh>
#include <ikarus/utils/polyfit.hh>

using Dune::TestSuite;

static auto polyFitTest1() {
  TestSuite t("polyFitTest1");
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(10, 0, 10);
  Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(10, 2, 20);

  auto [poly, normE] = Ikarus::utils::polyfit(x, y, 1);
  t.check(Dune::FloatCmp::eq(2.0, poly.coefficients()[0]));
  t.check(Dune::FloatCmp::eq(1.8, poly.coefficients()[1]));
  t.check(1e-14 > normE);
  return t;
}

static auto polyFitTest2() {
  TestSuite t("polyFitTest2");
  const double factor = 7.6;
  Eigen::VectorXd x   = Eigen::VectorXd::LinSpaced(10, 0, 10);
  Eigen::VectorXd y   = 7 * x.array().cwiseProduct(x.array()).matrix();
  for (int i = 0; i < y.size(); ++i) {
    y[i] += (1 - i / 10.0) * factor - (1 - i * i / 10.0) * factor + std::sin(i / 10.0);
  }

  auto [poly, normE] = Ikarus::utils::polyfit(x, y, 2);

  t.check(Dune::FloatCmp::eq(-0.0038062785674569739, poly.coefficients()[0]));
  t.check(Dune::FloatCmp::eq(-0.58760441700969401, poly.coefficients()[1]));
  t.check(Dune::FloatCmp::eq(7.6138682871655829, poly.coefficients()[2]));
  t.check(Dune::FloatCmp::eq(0.0082367593944499204, normE));
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(polyFitTest1());
  t.subTest(polyFitTest2());

  return t.exit();
}
