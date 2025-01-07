// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <experimental/type_traits>

#include <dune/common/test/testsuite.hh>
#include <dune/common/tuplevector.hh>
#if HAVE_DUNE_LOCALFEFUNCTIONS
  #include <dune/localfefunctions/manifolds/realTuple.hh>
  #include <dune/localfefunctions/manifolds/unitVector.hh>

#endif
#include <ikarus/solver/nonlinearsolver/trustregion.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/linearalgebrahelper.hh>
#include <ikarus/utils/nonlinearoperator.hh>

using namespace Ikarus;
using Dune::TestSuite;

static auto f(const Eigen::Vector<double, 1>& x) { return 0.5 * x[0] * x[0]; }
static auto df(const Eigen::Vector<double, 1>& x) {
  Eigen::Vector<double, 1> r;
  r[0] = x[0];
  return r;
}

static auto ddf(const Eigen::Vector<double, 1>&) {
  Eigen::SparseMatrix<double> A(1, 1);
  A.insert(0, 0);

  A.coeffRef(0, 0) = 1;

  return A;
}

static auto trustRegion1() {
  TestSuite t("trustRegion1");
  Eigen::Vector<double, 1> x;
  x << 2;

  auto fvLambda   = [](auto&& xL) { return f(xL); };
  auto dfvLambda  = [](auto&& xL) { return df(xL); };
  auto ddfvLambda = [](auto&& xL) { return ddf(xL); };
  NonLinearOperator nonLinOp(functions(fvLambda, dfvLambda, ddfvLambda), parameter(x));

  Eigen::Vector<double, 1> xExpected;
  xExpected << 0;

  auto tr = makeTrustRegion(nonLinOp);
  tr->setup({.verbosity = 1, .Delta0 = 1});
  const auto solverInfo = tr->solve();

  t.check(true == solverInfo.success);

  t.check(isApproxSame(x, xExpected, 1e-15));
  return t;
}

static constexpr double a_      = 1.0;
static constexpr double b_      = 100.0;
static constexpr double offset_ = -1;
static auto rosenbrock(const Eigen::Vector2d& x) {
  return Dune::power(a_ - x[0], 2) + b_ * Dune::power(x[1] - x[0] * x[0], 2) + offset_;
}
static auto rosenbrockdx(const Eigen::Vector2d& x) {
  Eigen::Vector2d r;
  r[0] = -2 * a_ + 2 * x[0] - 4 * b_ * (-x[0] * x[0] + x[1]) * x[0];
  r[1] = 2 * b_ * (-x[0] * x[0] + x[1]);
  return r;
}

static auto rosenbrockddx(const Eigen::Vector2d& x) {
  Eigen::SparseMatrix<double> A(2, 2);
  A.insert(0, 0);
  A.insert(0, 1);
  A.insert(1, 0);
  A.insert(1, 1);
  A.coeffRef(0, 0) = 2 + 8 * b_ * x[0] * x[0] - 4 * b_ * (-x[0] * x[0] + x[1]);
  A.coeffRef(0, 1) = A.coeffRef(1, 0) = -4 * b_ * x[0];
  A.coeffRef(1, 1)                    = 2 * b_;

  return A;
}

static auto trustRegion2() {
  TestSuite t("trustRegion2");
  Eigen::Vector2d x;
  x << 2, 3;

  auto fvLambda   = [](auto&& xL) { return rosenbrock(xL); };
  auto dfvLambda  = [](auto&& xL) { return rosenbrockdx(xL); };
  auto ddfvLambda = [](auto&& xL) { return rosenbrockddx(xL); };
  Ikarus::NonLinearOperator nonLinOp(functions(fvLambda, dfvLambda, ddfvLambda), parameter(x));
  const double eps   = 1e-10;
  const int maxIter_ = 30;
  Eigen::Vector2d xExpected;
  xExpected << a_, a_ * a_;

  Ikarus::TrustRegion tr(nonLinOp);
  tr.setup({.verbosity = 1, .maxIter = maxIter_, .grad_tol = eps, .Delta0 = 1});
  const auto solverInfo = tr.solve();

  t.check(true == solverInfo.success);
  t.check(25 == solverInfo.iterations);
  t.check(eps > solverInfo.residualNorm);
  t.check(isApproxSame(x, xExpected, eps));
  nonLinOp.update<0>();
  t.check(Dune::FloatCmp::eq(offset_, nonLinOp.value()));
  return t;
}

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
using namespace autodiff;

template <typename ScalarType>
ScalarType f3(const Eigen::Vector2<ScalarType>& x) {
  return -10 * x[0] * x[0] + 10 * x[1] * x[1] + 4 * sin(x[0] * x[1]) - 2 * x[0] + Dune::power(x[0], 4);
}
static Eigen::Vector2d df3(Eigen::Vector2<autodiff::dual>& x) {
  return autodiff::gradient(f3<autodiff::dual>, autodiff::wrt(x), autodiff::at(x));
}

static auto ddf3(Eigen::Vector2<autodiff::dual2nd>& x) {
  Eigen::SparseMatrix<double> A(2, 2);
  Eigen::Matrix2d h = autodiff::hessian(f3<autodiff::dual2nd>, autodiff::wrt(x), autodiff::at(x));
  A.insert(0, 0)    = h(0, 0);
  A.insert(0, 1)    = h(0, 1);
  A.insert(1, 0)    = h(1, 0);
  A.insert(1, 1)    = h(1, 1);

  return A;
}

static auto trustRegion3() {
  TestSuite t("trustRegion3");
  Eigen::Vector2d x(2);
  x << 0.7, -3.3;

  auto fvLambda  = [](auto&& xL) { return f3(xL); };
  auto dfvLambda = [](auto&& xL) {
    auto xR = xL.template cast<autodiff::dual>().eval();
    return df3(xR);
  };
  auto ddfvLambda = [](auto&& xL) {
    auto xR = xL.template cast<autodiff::dual2nd>().eval();
    return ddf3(xR);
  };
  Ikarus::NonLinearOperator nonLinOp(functions(fvLambda, dfvLambda, ddfvLambda), parameter(x));
  const double eps   = 1e-12;
  const int maxIter_ = 30;
  Eigen::Vector2d xExpected;
  xExpected << 2.3066301277034750861, -0.33230864873179355445;

  Ikarus::TrustRegion tr(nonLinOp);
  tr.setup({.verbosity = 1, .maxIter = maxIter_, .grad_tol = eps, .corr_tol = eps, .Delta0 = 1});
  const auto solverInfo = tr.solve();
  t.check(true == solverInfo.success);
  t.check(11 == solverInfo.iterations);
  t.check(eps > solverInfo.residualNorm);
  t.check(isApproxSame(x, xExpected, eps));
  nonLinOp.update<0>();
  t.check(Dune::FloatCmp::eq(-31.180733385187978, nonLinOp.value()));

  x << 0.7, -3.3;
  Ikarus::TrustRegion<decltype(nonLinOp), Ikarus::PreConditioner::IdentityPreconditioner> tr2(nonLinOp);
  tr2.setup({.verbosity = 1, .maxIter = maxIter_, .grad_tol = eps, .corr_tol = eps, .Delta0 = 1});
  const auto solverInfo2 = tr2.solve();
  t.check(true == solverInfo2.success);
  t.check(11 == solverInfo2.iterations);
  t.check(eps > solverInfo2.residualNorm);
  t.check(isApproxSame(x, xExpected, eps));
  nonLinOp.update<0>();
  t.check(Dune::FloatCmp::eq(-31.180733385187978, nonLinOp.value()));

  x << 0.7, -3.3;
  Ikarus::TrustRegion<decltype(nonLinOp), Ikarus::PreConditioner::DiagonalPreconditioner> tr3(nonLinOp);
  tr3.setup({.verbosity = 1, .maxIter = maxIter_, .grad_tol = eps, .corr_tol = eps, .Delta0 = 1});
  const auto solverInfo3 = tr3.solve();
  t.check(true == solverInfo3.success);
  t.check(8 == solverInfo3.iterations);
  t.check(eps > solverInfo3.residualNorm);
  t.check(isApproxSame(x, xExpected, eps));
  nonLinOp.update<0>();
  t.check(Dune::FloatCmp::eq(-31.180733385187978, nonLinOp.value()));
  return t;
}

template <typename ScalarType = double>
ScalarType f3R(const Dune::UnitVector<double, 2>& x,
               const Eigen::Vector<ScalarType, 2>& dx = Eigen::Vector<ScalarType, 2>::Zero()) {
  Eigen::Vector<ScalarType, 2> y = x.getValue();
  y += dx;
  return y[0] * y[0];
}
static Eigen::Vector<double, 1> df3R(const Dune::UnitVector<double, 2>& x) {
  Eigen::Vector<autodiff::dual, 2> xR = Eigen::Vector<autodiff::dual, 2>::Zero();
  auto dfvLambda                      = [&](auto&& xRL) { return f3R<autodiff::dual>(x, xRL); };
  autodiff::dual energy;
  Eigen::Vector<double, 2> g;
  autodiff::gradient(dfvLambda, autodiff::wrt(xR), autodiff::at(xR), energy, g);
  auto BLA = x.orthonormalFrame();
  return g.transpose() * BLA;
}

static auto ddf3R(const Dune::UnitVector<double, 2>& x) {
  Eigen::SparseMatrix<double> A(1, 1);
  Eigen::Vector<autodiff::dual2nd, 2> xR = Eigen::Vector<autodiff::dual2nd, 2>::Zero();
  auto dfvLambda                         = [&](auto&& xRL) { return f3R<autodiff::dual2nd>(x, xRL); };
  autodiff::dual2nd energy;
  Eigen::Vector2d g;
  Eigen::Matrix2d h;
  const auto y = x.getValue();
  autodiff::hessian(dfvLambda, autodiff::wrt(xR), autodiff::at(xR), energy, g, h);
  auto BLA       = x.orthonormalFrame();
  A.insert(0, 0) = BLA.transpose() * h * BLA - y.dot(g);

  return A;
}

static auto trustRegion4_RiemanianUnitSphere() {
  TestSuite t("trustRegion4_RiemanianUnitSphere");

  auto d = Dune::UnitVector<double, 2>(Eigen::Vector2d::UnitX() + 0.1 * Eigen::Vector2d::UnitY());
  d.update(Eigen::Vector<double, 1>::Ones());
  auto fvLambda = [](auto&& xL) { return f3R(xL); };

  auto dfvLambda  = [](auto&& xL) { return df3R(xL); };
  auto ddfvLambda = [](auto&& xL) { return ddf3R(xL); };

  Ikarus::NonLinearOperator nonLinOp(functions(fvLambda, dfvLambda, ddfvLambda), parameter(d));
  t.check(Dune::FloatCmp::eq(nonLinOp.value(), fvLambda(d))) << "Nonlinear operator and lambda have different value";

  t.check(isApproxSame(dfvLambda(d), nonLinOp.derivative(), 1e-15))
      << "Nonlinear operator derivative and lambda have different value";
  t.check(isApproxSame(ddfvLambda(d), nonLinOp.secondDerivative(), 1e-15))
      << "Nonlinear operator second derivative and lambda have different value";

  Ikarus::TrustRegion tr3(nonLinOp,
                          std::function([](Dune::UnitVector<double, 2>& x,
                                           const Dune::UnitVector<double, 2>::CorrectionType& d_) { x += d_; }));
  constexpr double tol = 1e-12;
  tr3.setup({.verbosity = 1, .maxIter = 1000, .grad_tol = tol, .corr_tol = tol, .Delta0 = 0.1});
  const auto solverInfo3 = tr3.solve();
  t.check(true == solverInfo3.success) << "Trust region was unsuccessful.";
  t.check(6 == solverInfo3.iterations) << "Trust region has not the expected numbers of iterations.";
  t.check(tol > solverInfo3.residualNorm) << "Trust region didn't reach the correct norm";
  nonLinOp.update<0>();
  t.check(1e-17 >= nonLinOp.value()) << "Trust region energy is not zero";
  t.check(isApproxSame(nonLinOp.firstParameter().getValue(), -Eigen::Vector2d::UnitY(), 1e-15))
      << "Trust region solution is wrong";
  return t;
}

#include <dune/istl/bvector.hh>
using DirectorVector     = Dune::BlockVector<Dune::UnitVector<double, 3>>;
using DisplacementVector = Dune::BlockVector<Dune::RealTuple<double, 3>>;
using MultiTypeVector    = Dune::TupleVector<DisplacementVector, DirectorVector>;

template <typename ScalarType = double>
ScalarType f3RBlocked(const MultiTypeVector& mT, const Eigen::VectorX<ScalarType>& dx) {
  using namespace Dune::Indices;
  auto& disp                          = mT[_0];
  auto& dir                           = mT[_1];
  Eigen::VectorX<ScalarType> dualDisp = Ikarus::viewAsFlatEigenVector(disp);
  Eigen::VectorX<ScalarType> dualDir  = Ikarus::viewAsFlatEigenVector(dir);
  dualDisp += dx.segment(0, dualDisp.size());
  dualDir += dx(Eigen::seq(dualDisp.size(), Eigen::indexing::last));
  ScalarType energy = 0;
  for (auto i = 0U; i < dir.size(); ++i) {
    energy += Dune::power(dualDir(i * 3), 2);
    energy += Dune::power(dualDir(i * 3 + 1), 2);
  }
  return energy + dualDisp.squaredNorm() + dualDisp.dot(dualDir);
}
static Eigen::VectorXd df3RBlocked(const MultiTypeVector& mT) {
  using namespace Dune::Indices;
  auto& disp = mT[_0];
  auto& dir  = mT[_1];
  auto dispE = Ikarus::viewAsFlatEigenVector(disp);
  auto dirE  = Ikarus::viewAsFlatEigenVector(dir);

  Eigen::VectorX<autodiff::dual> xR = Eigen::VectorXd::Zero(dispE.size() + dirE.size());
  auto dfvLambda                    = [&](auto&& xRL) { return f3RBlocked<autodiff::dual>(mT, xRL); };
  autodiff::dual energy;
  Eigen::VectorXd g;
  Eigen::VectorXd gRed(dispE.size() + dir.size() * dir[0].correctionSize);
  autodiff::gradient(dfvLambda, autodiff::wrt(xR), autodiff::at(xR), energy, g);

  gRed.segment(0, dispE.size()) = g.segment(0, dispE.size());
  for (auto i = 0U; i < dir.size(); ++i) {
    auto BLA                             = dir[i].orthonormalFrame();
    size_t indexStart                    = dispE.size() + i * dir[i].correctionSize;
    size_t indexStartE                   = dispE.size() + i * dir[i].valueSize;
    gRed.template segment<2>(indexStart) = BLA.transpose() * g.template segment<3>(indexStartE);
  }
  return gRed;
}

static auto ddf3RBlocked(const MultiTypeVector& mT) {
  Eigen::SparseMatrix<double> A;
  using namespace Dune::Indices;
  auto& disp = mT[_0];
  auto& dir  = mT[_1];
  auto dispE = Ikarus::viewAsFlatEigenVector(disp);
  auto dirE  = Ikarus::viewAsFlatEigenVector(dir);

  Eigen::VectorX<autodiff::dual2nd> xR = Eigen::VectorXd::Zero(dispE.size() + dirE.size());
  auto dfvLambda                       = [&](auto&& xRL) { return f3RBlocked<autodiff::dual2nd>(mT, xRL); };
  autodiff::dual2nd energy;
  Eigen::VectorXd g;
  Eigen::MatrixXd h;
  A.resize(dispE.size() + dir.size() * dir[0].correctionSize, dispE.size() + dir.size() * dir[0].correctionSize);
  autodiff::hessian(dfvLambda, autodiff::wrt(xR), autodiff::at(xR), energy, g, h);

  for (int i = 0; i < dispE.size(); ++i) {
    for (int j = 0; j < dispE.size(); ++j) {
      A.insert(i, j) = h(i, j);
    }
  }
  for (auto i = 0U; i < dir.size(); ++i) {
    Eigen::Index indexStartI  = dispE.size() + i * dir[0].correctionSize;
    Eigen::Index indexStartIE = dispE.size() + i * dir[0].valueSize;
    auto BLAI                 = dir[i].orthonormalFrame();
    for (auto j = 0U; j < dir.size(); ++j) {
      auto BLAJ = dir[j].orthonormalFrame();

      Eigen::Index indexStartJ  = dispE.size() + j * dir[0].correctionSize;
      Eigen::Index indexStartJE = dispE.size() + j * dir[0].valueSize;

      Eigen::Matrix2d blockIJ = BLAI.transpose() * h.block<3, 3>(indexStartIE, indexStartJE) * BLAJ;

      A.insert(indexStartI, indexStartJ)         = blockIJ(0, 0);
      A.insert(indexStartI, indexStartJ + 1)     = blockIJ(0, 1);
      A.insert(indexStartI + 1, indexStartJ)     = blockIJ(1, 0);
      A.insert(indexStartI + 1, indexStartJ + 1) = blockIJ(1, 1);
    }
    A.coeffRef(indexStartI, indexStartI) -= dir[i].getValue().dot(g.template segment<3>(indexStartIE));
    A.coeffRef(indexStartI + 1, indexStartI + 1) -= dir[i].getValue().dot(g.template segment<3>(indexStartIE));
  }

  for (auto i = 0U; i < disp.size(); ++i) {
    Eigen::Index indexStartI = i * disp[0].correctionSize;
    for (auto j = 0U; j < dir.size(); ++j) {
      auto BLAJ                                  = dir[j].orthonormalFrame();
      Eigen::Index indexStartJ                   = dispE.size() + j * dir[0].correctionSize;
      Eigen::Index indexStartJE                  = dispE.size() + j * dir[0].valueSize;
      Eigen::Matrix<double, 3, 2> blockIJ        = h.block<3, 3>(indexStartI, indexStartJE) * BLAJ;
      A.insert(indexStartI, indexStartJ)         = blockIJ(0, 0);
      A.insert(indexStartI + 1, indexStartJ)     = blockIJ(1, 0);
      A.insert(indexStartI + 2, indexStartJ)     = blockIJ(2, 0);
      A.insert(indexStartI, indexStartJ + 1)     = blockIJ(0, 1);
      A.insert(indexStartI + 1, indexStartJ + 1) = blockIJ(1, 1);
      A.insert(indexStartI + 2, indexStartJ + 1) = blockIJ(2, 1);
    }
  }

  for (auto i = 0U; i < dir.size(); ++i) {
    Eigen::Index indexStartI  = dispE.size() + i * dir[0].correctionSize;
    Eigen::Index indexStartIE = dispE.size() + i * dir[0].valueSize;
    auto BLAI                 = dir[i].orthonormalFrame();
    for (auto j = 0U; j < disp.size(); ++j) {
      Eigen::Index indexStartJ                   = j * disp[0].correctionSize;
      Eigen::Matrix<double, 2, 3> blockIJ        = BLAI.transpose() * h.block<3, 3>(indexStartIE, indexStartJ);
      A.insert(indexStartI, indexStartJ)         = blockIJ(0, 0);
      A.insert(indexStartI, indexStartJ + 1)     = blockIJ(0, 1);
      A.insert(indexStartI, indexStartJ + 2)     = blockIJ(0, 2);
      A.insert(indexStartI + 1, indexStartJ)     = blockIJ(1, 0);
      A.insert(indexStartI + 1, indexStartJ + 1) = blockIJ(1, 1);
      A.insert(indexStartI + 1, indexStartJ + 2) = blockIJ(1, 2);
    }
  }

  return A;
}

static auto trustRegion5_RiemanianUnitSphereAndDispBlocked() {
  TestSuite t("trustRegion5_RiemanianUnitSphereAndDispBlocked");
  using namespace Dune::Indices;
  using namespace Ikarus;
  DisplacementVector disp;
  disp.resize(2);
  disp[0].setValue(Eigen::Vector3d::UnitY());
  disp[1].setValue(Eigen::Vector3d::UnitZ());

  DirectorVector directors;
  directors.resize(2);
  directors[0].setValue(Eigen::Vector3d::UnitY() + Eigen::Vector3d::Ones());
  directors[1].setValue(Eigen::Vector3d::UnitZ() + Eigen::Vector3d::Ones());

  auto dispEigen     = Ikarus::viewAsFlatEigenVector(disp);
  auto directorEigen = Ikarus::viewAsFlatEigenVector(disp);

  Eigen::VectorXd zeroVec;
  zeroVec.setZero(dispEigen.size() + directorEigen.size());

  MultiTypeVector mT(disp, directors);

  auto fvLambda  = [&](auto&& xL) { return f3RBlocked(xL, zeroVec); };
  auto dfvLambda = [&](auto&& xL) { return df3RBlocked(xL); };

  auto gred = dfvLambda(mT);
  t.check(10 == gred.size());
  auto ddfvLambda = [](auto&& xL) { return ddf3RBlocked(xL); };
  auto h          = ddfvLambda(mT);
  t.check(10 == h.rows());
  t.check(10 == h.cols());

  Ikarus::NonLinearOperator nonLinOp(functions(fvLambda, dfvLambda, ddfvLambda), parameter(mT));
  t.check(Dune::FloatCmp::eq(nonLinOp.value(), fvLambda(mT)));

  t.check(isApproxSame(dfvLambda(mT), nonLinOp.derivative(), 1e-15));
  t.check(isApproxSame(ddfvLambda(mT), nonLinOp.secondDerivative(), 1e-15));

  TrustRegion tr3(nonLinOp);
  constexpr double tol = 1e-12;
  tr3.setup({.verbosity = 1, .maxIter = 1000, .grad_tol = tol, .corr_tol = tol, .Delta0 = 0.1});
  const auto solverInfo3 = tr3.solve();
  t.check(true == solverInfo3.success);
  t.check(9 == solverInfo3.iterations);
  t.check(tol > solverInfo3.residualNorm);
  nonLinOp.update<0>();
  t.check(Dune::FloatCmp::eq(-0.5, nonLinOp.value()));

  for (auto& director : mT[_1])
    t.check(isApproxSame(director.getValue(), Eigen::Vector3d::UnitZ(), 1e-15));
  for (auto& displacement : mT[_0])
    t.check(isApproxSame(displacement.getValue(), -0.5 * Eigen::Vector3d::UnitZ(), 1e-14));
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(trustRegion1());
  t.subTest(trustRegion2());
  t.subTest(trustRegion3());
  t.subTest(trustRegion4_RiemanianUnitSphere());
  t.subTest(trustRegion5_RiemanianUnitSphereAndDispBlocked());

  return t.exit();
}
