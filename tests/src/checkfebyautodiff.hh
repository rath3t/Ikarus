// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "testhelpers.hh"

#include <dune/common/float_cmp.hh>
#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/finiteelements/autodifffe.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/utils/basis.hh>

template <typename GridView, typename BasisHandler, typename Skills, typename AffordanceColl, typename VectorType>
auto checkFESByAutoDiffImpl(const GridView& gridView, const BasisHandler& basis, Skills&& skills,
                            AffordanceColl affordance, VectorType& d, const std::string& messageIfFailed = "",
                            double tol = 1e-10) {
  double lambda = 7.3;
  auto fe       = Ikarus::makeFE(basis, std::forward<Skills>(skills));
  using FE      = decltype(fe);

  auto req = typename FE::Requirement();
  req.insertGlobalSolution(d).insertParameter(lambda);
  Dune::TestSuite t("Check calculateScalarImpl() and calculateVectorImpl() by Automatic Differentiation" +
                    messageIfFailed);
  for (auto element : elements(gridView)) {
    auto localView = basis.flat().localView();
    localView.bind(element);
    auto nDOF = localView.size();

    fe.bind(element);

    const std::string feClassName = Dune::className(fe);

    using AutoDiffBasedFE = Ikarus::AutoDiffFE<decltype(fe), true>;
    AutoDiffBasedFE feAutoDiff(fe);

    Eigen::MatrixXd K, KAutoDiff;
    K.setZero(nDOF, nDOF);
    KAutoDiff.setZero(nDOF, nDOF);

    Eigen::VectorXd R, RAutoDiff;
    R.setZero(nDOF);
    RAutoDiff.setZero(nDOF);

    calculateMatrix(fe, req, affordance.matrixAffordance(), K);
    calculateMatrix(feAutoDiff, req, affordance.matrixAffordance(), KAutoDiff);

    calculateVector(fe, req, affordance.vectorAffordance(), R);
    calculateVector(feAutoDiff, req, affordance.vectorAffordance(), RAutoDiff);

    t.check(isApproxSame(K, KAutoDiff, tol),
            "Mismatch between the stiffness matrices obtained from explicit implementation and the one based on "
            "automatic differentiation for " +
                feClassName)
        << "K is \n"
        << K << "\nKAutoDiff is \n"
        << KAutoDiff << "\nThe difference is\n"
        << (K - KAutoDiff);

    t.check(isApproxSame(R, RAutoDiff, tol),
            "Mismatch between the residual vectors obtained from explicit implementation and the one based on "
            "automatic differentiation for " +
                feClassName)
        << "R is " << R.transpose() << "\nRAutoDiff is " << RAutoDiff.transpose() << "\nThe difference is "
        << (R - RAutoDiff).transpose();

    t.check(Dune::FloatCmp::eq(calculateScalar(fe, req, affordance.scalarAffordance()),
                               calculateScalar(feAutoDiff, req, affordance.scalarAffordance()), tol),
            "Mismatch between the energies obtained from explicit implementation and the one based on "
            "automatic differentiation for " +
                feClassName);
  }

  return t;
}

template <typename GridView, typename PreBasis, typename Skills, typename AffordanceColl>
auto checkFESByAutoDiff(const GridView& gridView, const PreBasis& pb, Skills&& skills, AffordanceColl affordance,
                        const std::string& testName = "", double tol = 1e-10) {
  Dune::TestSuite t("AutoDiff Test" + testName);
  auto basis = Ikarus::makeBasis(gridView, pb);
  Eigen::VectorXd d;
  d.setZero(basis.flat().dimension());
  t.subTest(checkFESByAutoDiffImpl(gridView, basis, skills, affordance, d, " Zero Displacements", tol));
  d.setRandom(basis.flat().dimension());
  d *= 0.2; // to avoid a negative determinant of deformation gradient
  t.subTest(checkFESByAutoDiffImpl(gridView, basis, skills, affordance, d, " Non-zero Displacements", tol));
  return t;
}
