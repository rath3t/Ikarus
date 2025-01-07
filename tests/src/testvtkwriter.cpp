// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "dummyproblem.hh"
#include "testhelpers.hh"

#include <memory>
#include <type_traits>

#include "dune/common/indices.hh"
#include "dune/functions/functionspacebases/subspacebasis.hh"
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/common/tuplevector.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/vtk/datacollectors/discontinuousdatacollector.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/vtkreader.hh>

#include <Eigen/Core>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/feresulttypes.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mechanics/loads/volume.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/io/resultfunction.hh>
#include <ikarus/io/vtkwriter.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

auto vtkWriterTest() {
  TestSuite t("Test ResultFunction");
  std::string fileName = "ResultFunctionTest";

  using Grid     = Dune::UGGrid<2>;
  using GridView = Grid::LeafGridView;

  DummyProblem<Grid> testCase{};

  auto& gridView       = testCase.gridView();
  auto sparseAssembler = testCase.sparseAssembler();
  auto& req            = testCase.requirement();
  auto& basis          = testCase.basis();
  auto& D_Glob         = req.globalSolution();

  // Tests
  Dune::Vtk::DiscontinuousDataCollector dc{gridView};

  auto writer = Ikarus::Vtk::Writer(sparseAssembler);

  using Ikarus::Vtk::DataTag::asCellAndPointData;
  using Ikarus::Vtk::DataTag::asCellData;
  using Ikarus::Vtk::DataTag::asPointData;

  writer.setDatatype(Dune::Vtk::DataTypes::FLOAT64);
  writer.setFormat(Dune::Vtk::FormatTypes::ASCII);

  writer.addResult<Ikarus::ResultTypes::linearStress>(); // Defaults to pointData
  writer.addResultFunction(Ikarus::makeResultFunction<Ikarus::ResultTypes::linearStress>(sparseAssembler), asCellData);

  writer.addInterpolation(D_Glob, basis.flat(), "displacement", asPointData);

  auto subspaceBasis = Dune::Functions::subspaceBasis(basis.flat(), Dune::index_constant<0>());
  writer.addInterpolation(D_Glob, subspaceBasis, "displacement_u", asCellAndPointData);

  writer.addPointData(
      Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis.flat(), D_Glob),
      Dune::Vtk::FieldInfo("displacements_gf", 2, Dune::Vtk::RangeTypes::VECTOR));

  auto vtkFileName = writer.write(fileName);

  Dune::GridFactory<Grid> factory;
  Dune::Vtk::VtkReader reader{factory};

  reader.read(vtkFileName);
  auto gridRead     = reader.createGrid();
  auto gridViewRead = gridRead->leafGridView();

  t.check(gridView.size(0) == gridViewRead.size(0));

  auto stressData = reader.getPointData("linearStress");
  t.check(stressData.numComponents() == 3)
      << testLocation() << "Num components should be 3, but is " << stressData.numComponents();
  t.check(stressData.dataType() == Dune::Vtk::DataTypes::FLOAT64)
      << testLocation() << std::source_location::current().line() << "Precision should be float64";

  auto stressDataCell = reader.getCellData("linearStress");
  t.check(stressDataCell.numComponents() == 3)
      << testLocation() << "Num components should be 3, but is " << stressDataCell.numComponents();
  t.check(stressDataCell.dataType() == Dune::Vtk::DataTypes::FLOAT64)
      << testLocation() << "Precision is should be float64";

  // Writing out as a vector does always seem to yield 3 (or potentially more) entries
  auto displacementData = reader.getPointData("displacement");
  t.check(displacementData.numComponents() == 3)
      << testLocation() << "Num components should be 3, but is " << displacementData.numComponents();

  auto displacementDataGF = reader.getPointData("displacements_gf");
  t.check(displacementDataGF.numComponents() == 3)
      << testLocation() << "Num components should be 3, but is " << displacementDataGF.numComponents();

  auto displacementDataU = reader.getPointData("displacement_u");
  t.check(displacementDataU.numComponents() == 1)
      << testLocation() << "Num components should be 1, but is " << displacementDataU.numComponents();

  auto writer2 = Ikarus::Vtk::Writer(sparseAssembler, Dune::Vtk::DiscontinuousDataCollector<GridView>{gridView});
  writer2.addAllResults(asPointData);
  writer2.addAllResults(asCellData);
  auto vtkFileName2 = writer2.write(fileName + "_2");

  reader.read(vtkFileName2);
  t.checkNoThrow([&]() { reader.getPointData("linearStress"); });
  t.checkNoThrow([&]() { reader.getCellData("linearStress"); });

  return t;
}

auto testUnstructuredInstantiaionAndDeduction() {
  TestSuite t;

  using Grid     = Dune::UGGrid<2>;
  using GridView = Grid::LeafGridView;

  DummyProblem<Grid> testCase{};

  auto& gridView       = testCase.gridView();
  auto sparseAssembler = testCase.sparseAssembler();

  static_assert(not Ikarus::Vtk::IsStructured<Grid>::value);

  Ikarus::Vtk::Writer writer(sparseAssembler);

  using DC   = typename decltype(writer)::DataCollector;
  using VTKW = typename decltype(writer)::VTKWriter;

  static_assert(std::is_same_v<DC, Dune::Vtk::ContinuousDataCollector<GridView>>);
  static_assert(std::is_same_v<VTKW, Dune::Vtk::UnstructuredGridWriter<GridView, DC>>);

  // Second writer with a DiscontinousWriter
  auto discontinousDC = Dune::Vtk::DiscontinuousDataCollector<GridView>(gridView);
  Ikarus::Vtk::Writer writer2(sparseAssembler, discontinousDC);

  using DC2   = typename decltype(writer2)::DataCollector;
  using VTKW2 = typename decltype(writer2)::VTKWriter;

  static_assert(std::is_same_v<DC2, Dune::Vtk::DiscontinuousDataCollector<GridView>>);
  static_assert(std::is_same_v<VTKW2, Dune::Vtk::UnstructuredGridWriter<GridView, decltype(discontinousDC)>>);

  // Use Args
  Ikarus::Vtk::Writer writer4(sparseAssembler, Dune::Vtk::FormatTypes::BINARY, Dune::Vtk::DataTypes::FLOAT32);
  Ikarus::Vtk::Writer writer5(sparseAssembler, discontinousDC, Dune::Vtk::FormatTypes::BINARY,
                              Dune::Vtk::DataTypes::FLOAT32);
  Ikarus::Vtk::Writer writer6(sparseAssembler, std::move(discontinousDC), Dune::Vtk::FormatTypes::BINARY,
                              Dune::Vtk::DataTypes::FLOAT32);

  return t;
}

auto testStructuredInstantiaionAndDeduction() {
  TestSuite t;

  using Grid     = Dune::YaspGrid<2>;
  using GridView = Grid::LeafGridView;

  DummyProblem<Grid, true> testCase{};

  auto& gridView       = testCase.gridView();
  auto sparseAssembler = testCase.sparseAssembler();

  static_assert(Ikarus::Vtk::IsStructured<Grid>::value);

  Ikarus::Vtk::Writer writer(sparseAssembler);

  using DC   = typename decltype(writer)::DataCollector;
  using VTKW = typename decltype(writer)::VTKWriter;

  static_assert(std::is_same_v<DC, Dune::Vtk::YaspDataCollector<GridView>>);
  static_assert(std::is_same_v<VTKW, Dune::Vtk::RectilinearGridWriter<GridView, DC>>);

  // Second writer with a StructuredDataCollector
  auto structuredDC = Dune::Vtk::StructuredDataCollector<GridView>(gridView);
  Ikarus::Vtk::Writer writer2(sparseAssembler, structuredDC);

  using DC2   = typename decltype(writer2)::DataCollector;
  using VTKW2 = typename decltype(writer2)::VTKWriter;

  static_assert(std::is_same_v<DC2, Dune::Vtk::StructuredDataCollector<GridView>>);
  static_assert(std::is_same_v<VTKW2, Dune::Vtk::RectilinearGridWriter<GridView, decltype(structuredDC)>>);

  // Use Args
  Ikarus::Vtk::Writer writer4(sparseAssembler, Dune::Vtk::FormatTypes::BINARY, Dune::Vtk::DataTypes::FLOAT32);
  Ikarus::Vtk::Writer writer5(sparseAssembler, structuredDC, Dune::Vtk::FormatTypes::BINARY,
                              Dune::Vtk::DataTypes::FLOAT32);
  Ikarus::Vtk::Writer writer6(sparseAssembler, std::move(structuredDC), Dune::Vtk::FormatTypes::BINARY,
                              Dune::Vtk::DataTypes::FLOAT32);

  return t;
}

int main(const int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testStructuredInstantiaionAndDeduction());
  t.subTest(testUnstructuredInstantiaionAndDeduction());
  t.subTest(vtkWriterTest());

  return t.exit();
}