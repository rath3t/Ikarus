// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <ikarus/utils/autodiffhelper.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

template <typename Scalar>
Eigen::Vector<Scalar, 2> f(const Eigen::Vector<Scalar, 3>& x) {
  return (x.array() * (x.array().sin())).template segment<2>(0);
}

static auto hessianN() {
  TestSuite t("hessianN");
  Eigen::Vector3d xd;
  xd << 1.0, 2.0, 3.0;
  Eigen::Vector3dual2nd x = xd;
  Eigen::Vector2dual2nd u;
  std::array<Eigen::Vector<double, 3>, 2> g;
  std::array<Eigen::Matrix<double, 3, 3>, 2> h;
  Ikarus::utils::hessianN(f<autodiff::dual2nd>, wrt(x), at(x), u, g, h);

  for (Eigen::Index i = 0; i < 2; ++i) {
    Eigen::Vector3d gExpected;
    Eigen::Matrix3d hExpected;
    gExpected.setZero();
    hExpected.setZero();
    gExpected[i]    = sin(xd[i]) + cos(xd[i]) * xd[i];
    hExpected(i, i) = 2 * cos(xd[i]) - xd[i] * sin(xd[i]);
    t.check(isApproxSame(g[static_cast<size_t>(i)], gExpected, 1e-14));
    t.check(isApproxSame(h[static_cast<size_t>(i)], hExpected, 1e-14));
  }
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(hessianN());

  return t.exit();
}
