// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#include <config.h>

#include "dummyproblem.hh"
#include "testhelpers.hh"

#include <dune/common/float_cmp.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mechanics/loads/volume.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/solver/eigenvaluesolver/generalizedeigensolver.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

template <Ikarus::Concepts::EigenValueSolver SOL1, Ikarus::Concepts::EigenValueSolver SOL2>
auto testEigenValues(const SOL1& solver1, const SOL2& solver2, std::optional<Eigen::Index> nev_ = std::nullopt) {
  TestSuite t("Test Eigenvalues" + Dune::className<SOL1>() + ", " + Dune::className<SOL2>());
  auto eigenvalues1 = solver1.eigenvalues();
  auto eigenvalues2 = solver2.eigenvalues();

  if (not nev_.has_value()) {
    t.check(isApproxSame(eigenvalues1, eigenvalues2, 1e-10)) << testLocation() << "\n"
                                                             << eigenvalues1 << "\n\n"
                                                             << eigenvalues2;
  } else {
    t.check(isApproxSame(eigenvalues1.head(*nev_).eval(), eigenvalues2.head(*nev_).eval(), 1e-10))
        << testLocation() << "\n"
        << eigenvalues1.head(*nev_).eval() << "\n\n"
        << eigenvalues2.head(*nev_).eval();
  }

  return t;
}

// as1 corresponds to solver1, which is tested against solver2
template <Ikarus::Concepts::EigenValueSolver SOL1, Ikarus::Concepts::EigenValueSolver SOL2,
          Ikarus::Concepts::FlatAssembler AS1>
auto testEigenVectors(const SOL1& solver1, const SOL2& solver2, std::shared_ptr<AS1> as1) {
  TestSuite t("Test Eigenvalues " + Dune::className<SOL1>() + ", " + Dune::className<SOL2>());

  Eigen::MatrixXd K = as1->matrix();

  Eigen::MatrixXd evecs = solver1.eigenvectors();
  Eigen::VectorXd evals = solver2.eigenvalues();

  Eigen::VectorXd T = (evecs.transpose() * K * evecs).diagonal();

  t.check(isApproxSame(T, evals, 1e-10)) << testLocation() << "\n" << T << "\n\n" << evals;

  return t;
}

auto testEigenSolvers() {
  TestSuite t("EigenValueSolver Test");
  using Grid     = Dune::YaspGrid<2>;
  using GridView = Grid::LeafGridView;

  DummyProblem<Grid, true> testCase({5, 5});
  auto& req             = testCase.requirement();
  auto& fes             = testCase.finiteElements();
  auto& dirichletValues = testCase.dirichletValues();

  auto assM = Ikarus::makeSparseFlatAssembler(fes, dirichletValues);
  auto assK = Ikarus::makeSparseFlatAssembler(fes, dirichletValues);

  auto assMD = Ikarus::makeDenseFlatAssembler(fes, dirichletValues);
  auto assKD = Ikarus::makeDenseFlatAssembler(fes, dirichletValues);

  assM->bind(req, Ikarus::AffordanceCollections::modalAnalysis, Ikarus::DBCOption::Reduced);
  assK->bind(req, Ikarus::AffordanceCollections::elastoStatics, Ikarus::DBCOption::Reduced);

  assMD->bind(req, Ikarus::AffordanceCollections::modalAnalysis, Ikarus::DBCOption::Reduced);
  assKD->bind(req, Ikarus::AffordanceCollections::elastoStatics, Ikarus::DBCOption::Reduced);

  int nev = 10; // number of requested eigenvalues
  using Ikarus::EigenValueSolverType;

  auto partialSolver = Ikarus::PartialGeneralizedSymEigenSolver(assK, assM, nev);
  t.checkThrow([&]() { partialSolver.eigenvalues(); }) << testLocation();
  t.check(partialSolver.compute()) << testLocation();

  auto partialSolverDense = Ikarus::PartialGeneralizedSymEigenSolver(assKD, assMD, nev);
  t.check(partialSolverDense.compute()) << testLocation();

  auto solver1 = Ikarus::makeGeneralizedSymEigenSolver<EigenValueSolverType::Eigen>(assKD, assMD);
  t.check(solver1.compute()) << testLocation();

  auto solver2 = Ikarus::makeGeneralizedSymEigenSolver<EigenValueSolverType::Spectra>(assK, assM);
  t.check(solver2.compute()) << testLocation();

  auto solver3 = Ikarus::makeGeneralizedSymEigenSolver<EigenValueSolverType::Spectra>(assKD, assMD);
  t.check(solver3.compute()) << testLocation();

  t.subTest(testEigenVectors(solver2, solver3, assK));
  t.subTest(testEigenVectors(solver3, solver2, assKD));
  t.subTest(testEigenVectors(solver1, solver2, assKD));
  t.subTest(testEigenVectors(solver2, solver1, assK));

  t.subTest(testEigenValues(partialSolver, solver1, nev));
  t.subTest(testEigenValues(partialSolver, solver2, nev));
  t.subTest(testEigenValues(partialSolver, solver3, nev));
  t.subTest(testEigenValues(solver1, solver2));
  t.subTest(testEigenValues(solver1, solver3));
  t.subTest(testEigenValues(solver2, solver3));
  t.subTest(testEigenValues(partialSolver, partialSolverDense));

  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(partialSolver)>);
  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver1)>);
  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver2)>);
  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver3)>);

  auto partialIdentitySolver1 = Ikarus::makePartialIdentitySymEigenSolver(assK, 10);
  auto partialIdentitySolver2 = Ikarus::makePartialIdentitySymEigenSolver(assKD, 10);

  t.check(partialIdentitySolver1.compute()) << testLocation();
  t.check(partialIdentitySolver2.compute()) << testLocation();

  t.subTest(testEigenValues(partialIdentitySolver1, partialIdentitySolver2));

  auto solver1Identity = Ikarus::makeGeneralizedSymEigenSolver<EigenValueSolverType::Eigen>(assKD, assMD);
  t.check(solver1Identity.compute()) << testLocation();

  auto solver2Identity = Ikarus::makeGeneralizedSymEigenSolver<EigenValueSolverType::Spectra>(assK, assM);
  t.check(solver2Identity.compute()) << testLocation();

  auto solver3Identity = Ikarus::makeGeneralizedSymEigenSolver<EigenValueSolverType::Spectra>(assKD, assMD);
  t.check(solver3Identity.compute()) << testLocation();

  t.subTest(testEigenVectors(solver2Identity, solver3Identity, assK));
  t.subTest(testEigenVectors(solver3Identity, solver2Identity, assKD));
  t.subTest(testEigenVectors(solver1, solver2Identity, assKD));
  t.subTest(testEigenVectors(solver2Identity, solver1Identity, assK));

  t.subTest(testEigenValues(solver1Identity, solver2Identity));
  t.subTest(testEigenValues(solver1Identity, solver3Identity));
  t.subTest(testEigenValues(solver2Identity, solver3Identity));

  return t;
}

auto checkThrows() {
  TestSuite t;
  using Ikarus::EigenValueSolverType;

  auto testMatrix1 = Eigen::MatrixX<double>::Random(10, 10).eval();
  auto testMatrix2 = Eigen::MatrixX<double>::Random(9, 9).eval();
  auto testMatrix3 = Eigen::MatrixX<double>::Random(9, 10).eval();

  auto sparseTestMatrix1 = Eigen::SparseMatrix<double>(testMatrix1.sparseView());
  auto sparseTestMatrix2 = Eigen::SparseMatrix<double>(testMatrix2.sparseView());
  auto sparseTestMatrix3 = Eigen::SparseMatrix<double>(testMatrix3.sparseView());

  t.checkThrow([&]() {
    Ikarus::GeneralizedSymEigenSolver<EigenValueSolverType::Spectra, Eigen::SparseMatrix<double>>(sparseTestMatrix1,
                                                                                                  sparseTestMatrix2);
  }) << testLocation();
  t.checkThrow([&]() {
    Ikarus::GeneralizedSymEigenSolver<EigenValueSolverType::Spectra, Eigen::SparseMatrix<double>>(sparseTestMatrix2,
                                                                                                  sparseTestMatrix1);
  }) << testLocation();
  t.checkThrow([&]() {
    Ikarus::GeneralizedSymEigenSolver<EigenValueSolverType::Spectra, Eigen::SparseMatrix<double>>(sparseTestMatrix3,
                                                                                                  sparseTestMatrix1);
  }) << testLocation();

  t.checkThrow([&]() {
    Ikarus::GeneralizedSymEigenSolver<EigenValueSolverType::Eigen, Eigen::MatrixX<double>>(testMatrix1, testMatrix2);
  }) << testLocation();
  t.checkThrow([&]() {
    Ikarus::GeneralizedSymEigenSolver<EigenValueSolverType::Spectra, Eigen::MatrixX<double>>(testMatrix2, testMatrix1);
  }) << testLocation();
  t.checkThrow([&]() {
    Ikarus::GeneralizedSymEigenSolver<EigenValueSolverType::Spectra, Eigen::MatrixX<double>>(testMatrix3, testMatrix1);
  }) << testLocation();

  t.checkThrow([&]() { Ikarus::PartialGeneralizedSymEigenSolver<Eigen::MatrixX<double>>(testMatrix1, testMatrix2, 7); })
      << testLocation();
  t.checkThrow([&]() { Ikarus::PartialGeneralizedSymEigenSolver<Eigen::MatrixX<double>>(testMatrix2, testMatrix1, 7); })
      << testLocation();
  t.checkThrow([&]() { Ikarus::PartialGeneralizedSymEigenSolver<Eigen::MatrixX<double>>(testMatrix3, testMatrix1, 7); })
      << testLocation();

  t.checkThrow([&]() {
    Ikarus::PartialGeneralizedSymEigenSolver<Eigen::MatrixX<double>>(testMatrix1, testMatrix1, 15);
  }) << testLocation();

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testEigenSolvers());
  t.subTest(checkThrows());

  return t.exit();
}
