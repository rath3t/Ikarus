// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#include <config.h>

#include "tests/src/testhelpers.hh"

#include <vector>

#include <dune/common/float_cmp.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/fufem/dunepython.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/localfefunctions/eigenDuneTransformations.hh>
#include <dune/vtk/datacollectors/discontinuousdatacollector.hh>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mechanics/loads/volume.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/io/resultfunction.hh>
#include <ikarus/io/vtkwriter.hh>
#include <ikarus/solver/eigenvaluesolver/generaleigensolver.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/dynamics/dynamics.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

static auto dynamicsTest() {
  TestSuite t("DynamicsTest");
  using Grid = Dune::YaspGrid<2>;

  const double Lx                         = 4.0;
  const double Ly                         = 4.0;
  Dune::FieldVector<double, 2> bbox       = {Lx, Ly};
  std::array<int, 2> elementsPerDirection = {4, 4};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);

  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>()));

  auto vL      = []([[maybe_unused]] auto& globalCoord, auto& lamb) { return Eigen::Vector2d{0, -1}; };
  auto linMat  = Ikarus::LinearElasticity(Ikarus::toLamesFirstParameterAndShearModulus({.emodul = 100, .nu = 0.2}));
  auto skills_ = Ikarus::skills(Ikarus::linearElastic(Ikarus::planeStress(linMat)), Ikarus::volumeLoad<2>(vL));

  using FEType = decltype(Ikarus::makeFE(basis, skills_));
  std::vector<FEType> fes;
  for (auto&& element : elements(gridView)) {
    fes.emplace_back(Ikarus::makeFE(basis, skills_));
    fes.back().bind(element);
  }

  auto& fe = fes[0];
  auto req = FEType::Requirement();

  typename FEType::Requirement::SolutionVectorType d;
  d.setZero(basis.flat().size());
  req.insertGlobalSolution(d);

  Eigen::MatrixX<double> m(8, 8);
  m.setZero();
  fe.calculateMatrixImpl<double>(req, Ikarus::MatrixAffordance::mass, m);

  auto dirichletValues = Ikarus::DirichletValues(basis.flat());
  dirichletValues.fixBoundaryDOFs([](auto& dirichFlags, auto&& indexGlobal) { dirichFlags[indexGlobal] = true; });

  auto assM = Ikarus::makeSparseFlatAssembler(fes, dirichletValues);
  auto assK = Ikarus::makeSparseFlatAssembler(fes, dirichletValues);

  auto assMD = Ikarus::makeDenseFlatAssembler(fes, dirichletValues);
  auto assKD = Ikarus::makeDenseFlatAssembler(fes, dirichletValues);

  assM->bind(req, Ikarus::AffordanceCollections::dynamics, Ikarus::DBCOption::Reduced);
  assK->bind(req, Ikarus::AffordanceCollections::elastoStatics, Ikarus::DBCOption::Reduced);

  assMD->bind(req, Ikarus::AffordanceCollections::dynamics, Ikarus::DBCOption::Reduced);
  assKD->bind(req, Ikarus::AffordanceCollections::elastoStatics, Ikarus::DBCOption::Reduced);

  auto assMLumped = Ikarus::Dynamics::makeLumpedFlatAssembler(assM);

  auto mat1 = assMLumped->matrix();
  auto mat2 = assMLumped->matrix();

  int nev = 10; // number of requested eigenvalues
  using Ikarus::EigenSolverTypeTag;

  auto solver  = Ikarus::PartialGeneralSymEigenSolver(assK, assM, nev);
  bool success = solver.compute();
  std::cout << "Success: " << std::boolalpha << success << std::endl;

  auto eigenvectors = solver.eigenvectors();
  auto eigenvalues  = solver.eigenvalues();

  std::cout << "Angular frequencies found (n=" << nev << "):\n" << eigenvalues.array().sqrt() << std::endl;

  auto solver2 = Ikarus::PartialGeneralSymEigenSolver(assK, assMLumped, nev);
  solver2.compute();

  auto eigenvalues2 = solver2.eigenvalues();
  std::cout << "Angular frequencies found (n=" << nev << "):\n" << eigenvalues2.array().sqrt() << std::endl;

  Ikarus::Dynamics::writeEigenformsToVTK(solver2, assM, "eigenforms");
  Ikarus::Dynamics::writeEigenformsToPVD(solver, assM, "eigenforms_full");
  // Ikarus::Dynamics::writeEigenformsToPVD(solver2, assM, "eigenforms_lumped");

  auto solver3 = Ikarus::makeGeneralSymEigenSolver<EigenSolverTypeTag::Eigen>(assKD, assMD);
  solver3.compute();
  auto eigenvalues3 = solver3.eigenvalues();

  auto solver4 = Ikarus::makeGeneralSymEigenSolver<EigenSolverTypeTag::Spectra>(assK, assM);
  solver4.compute();
  auto eigenvalues4 = solver4.eigenvalues();

  auto solver5 = Ikarus::makeGeneralSymEigenSolver<EigenSolverTypeTag::Spectra>(assKD, assMD);
  solver5.compute();
  auto eigenvalues5 = solver5.eigenvalues();

  t.check(isApproxSame(eigenvalues3, eigenvalues4, 1e-10));
  t.check(isApproxSame(eigenvalues5, eigenvalues4, 1e-10));
  t.check(isApproxSame(eigenvalues4.head(nev).eval(), eigenvalues, 1e-10));

  Ikarus::Dynamics::writeEigenformsToVTK(solver4, assM, "eigenforms_vtk");
  Ikarus::Dynamics::writeEigenformsToPVD(solver4, assM, "eigenforms_pvd");

  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver)>);
  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver2)>);
  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver3)>);
  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver4)>);
  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver5)>);

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(dynamicsTest());
  return t.exit();
}