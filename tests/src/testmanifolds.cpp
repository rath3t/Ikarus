// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <array>

#include <dune/common/test/testsuite.hh>
#if HAVE_DUNE_LOCALFEFUNCTIONS
  #include <dune/localfefunctions/manifolds/realTuple.hh>
  #include <dune/localfefunctions/manifolds/unitVector.hh>
#endif
#include <Eigen/Core>

#include <ikarus/utils/init.hh>

using Dune::TestSuite;

static constexpr double tol = 1e-15;

static auto testUnitVector() {
  TestSuite t("testUnitVector");
  using namespace Dune;
  UnitVector<double, 3> a{UnitVector<double, 3>::CoordinateType::UnitZ()};
  a.update(Eigen::Vector<double, 2>::UnitX());
  const auto aExpected = Eigen::Vector<double, 3>(1.0 / sqrt(2), 0.0, 1.0 / sqrt(2));
  t.check(isApproxSame(a.getValue(), aExpected, tol));

  a = update(a, Eigen::Vector<double, 2>::UnitY());
  t.check(isApproxSame(a.getValue(), Eigen::Vector<double, 3>(0.5, 1.0 / sqrt(2), 0.5), tol));

  const auto d = update(a, Eigen::Vector<double, 2>::UnitY());

  t.check(isApproxSame(
      d.getValue(), Eigen::Vector<double, 3>(0.18688672392660707344, 0.97140452079103167815, -0.14644660940672621363),
      tol));

  UnitVector<double, 3> b{a};

  t.check(isApproxSame(b.getValue(), Eigen::Vector<double, 3>(0.5, 1.0 / sqrt(2), 0.5), tol));

  UnitVector<double, 3> c{UnitVector<double, 3>{Eigen::Vector<double, 3>::UnitZ() * 2.0}}; // move constructor test
  t.check(isApproxSame(c.getValue(), Eigen::Vector<double, 3>(0.0, 0.0, 1.0), tol));

  c.setValue(Eigen::Vector<double, 3>(13.0, -5.0, 1.0));
  t.check(isApproxSame(c.getValue(), Eigen::Vector<double, 3>(13.0, -5.0, 1.0).normalized(), tol));

  b = a;
  t.check(a == b);
  const auto testVec = Eigen::Vector<double, 3>(127.0, -5.0, 1.0);
  b.setValue(testVec);

  t.check(isApproxSame(b.getValue(), testVec.normalized(), tol));

  auto e{std::move(a)};

  e.setValue(Eigen::Vector<double, 3>(0.0, 0.0, -1.0));

  e = update(e, Eigen::Vector<double, 2>::UnitY());

  const auto eExpected = Eigen::Vector<double, 3>(0, 1.0 / sqrt(2), -1.0 / sqrt(2));
  t.check(isApproxSame(e.getValue(), eExpected, tol));
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testUnitVector());

  return t.exit();
}
