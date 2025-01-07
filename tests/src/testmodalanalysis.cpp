// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#include <config.h>

#include "tests/src/dummyproblem.hh"
#include "tests/src/testhelpers.hh"

#include <dune/common/float_cmp.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/modalanalysis/modalanalysis.hh>

using Dune::TestSuite;

static auto modalAnalysisTest() {
  TestSuite t("ModalAnalysis Test");
  using Grid = Dune::YaspGrid<2>;

  DummyProblem<Grid, true> testCase({5, 5});
  auto& req             = testCase.requirement();
  auto& fes             = testCase.finiteElements();
  auto& dirichletValues = testCase.dirichletValues();
  auto mA               = Ikarus::Dynamics::ModalAnalysis(std::move(fes), dirichletValues);
  mA.compute();

  auto frequencies = mA.angularFrequencies();

  mA.bindLumpingScheme<Ikarus::Dynamics::LumpingSchemes::RowSumLumping>();
  mA.compute();

  auto frequenciesLumped = mA.angularFrequencies();
  t.check(frequencies.sum() > frequenciesLumped.sum()) << testLocation();

  mA.unBindLumpingSchemes();
  mA.compute();
  auto frequencies2 = mA.angularFrequencies();
  t.check(isApproxSame(frequencies, frequencies2, 1e-14)) << testLocation();

  mA.writeEigenModes("eigenforms", 20);

  auto massAssembler = Ikarus::makeSparseFlatAssembler(fes, dirichletValues);
  massAssembler->bind(req, Ikarus::AffordanceCollections::modalAnalysis, Ikarus::DBCOption::Reduced);
  auto lumpedAssembler = Ikarus::Dynamics::makeLumpedFlatAssembler(massAssembler);

  auto M       = massAssembler->matrix();
  auto MLumped = lumpedAssembler->matrix();

  t.check(MLumped.nonZeros() == massAssembler->reducedSize()) << testLocation();
  t.check(MLumped.coeff(0, 0) == M.row(0).sum()) << testLocation();

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(modalAnalysisTest());
  return t.exit();
}
