// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "dummyproblem.hh"
#include "testhelpers.hh"

#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/vtk/datacollectors/discontinuousdatacollector.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/vtkreader.hh>

#include "ikarus/finiteelements/feresulttypes.hh"
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/io/resultfunction.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

auto expectedValues(Dune::Vtk::DataTypes precision, bool native) -> std::vector<double> {
  if (precision == Dune::Vtk::DataTypes::FLOAT32)
    return {-0.612563, -3.06282,   -0.178688,  -0.605305,  -3.02653,   1.40567e-16, -0.165843,   -2.97347,
            -0.164172, -0.158585,  -2.93718,   0.0145161,  -0.605305,  -3.02653,    1.40567e-16, -0.612563,
            -3.06282,  0.178688,   -0.158585,  -2.93718,   -0.0145161, -0.165843,   -2.97347,    0.164172,
            0.250011,  -0.894201,  0.122203,   0.218461,   -1.05195,   0.0145161,   -0.0192071,  -0.948045,
            0.0591021, -0.0507578, -1.1058,    -0.0485852, 0.218461,   -1.05195,    -0.0145161,  0.250011,
            -0.894201, -0.122203,  -0.0507578, -1.1058,    0.0485852,  -0.0192071,  -0.948045,   -0.0591021};

  if (precision == Dune::Vtk::DataTypes::FLOAT64 and native)
    return {-0.612563441928246,   -3.06281720964123,  -0.178688142264123,  -0.605305372298166,  -3.02652686149083,
            1.40566543035979e-16, -0.165843086267937, -2.97347313850917,   -0.164172003003964,  -0.158585016637857,
            -2.93718279035877,    0.0145161392601598, -0.605305372298166,  -3.02652686149083,   1.40566543035979e-16,
            -0.612563441928246,   -3.06281720964123,  0.178688142264124,   -0.158585016637857,  -2.93718279035877,
            -0.0145161392601598,  -0.165843086267937, -2.97347313850917,   0.164172003003964,   0.25001124825652,
            -0.894201465886881,   0.122203474568825,  0.218460568142139,   -1.05195486645879,   0.0145161392601599,
            -0.0192070900151417,  -0.948045133541213, 0.0591021143400626,  -0.0507577701295228, -1.10579853411312,
            -0.0485852209686023,  0.218460568142139,  -1.05195486645879,   -0.0145161392601597, 0.25001124825652,
            -0.894201465886881,   -0.122203474568824, -0.0507577701295225, -1.10579853411312,   0.0485852209686022,
            -0.0192070900151415,  -0.948045133541214, -0.0591021143400625};

  return {
      -0.612563441928246,    -3.06281720964123,    -0.1786881422641234,  -0.6053053722981661,  -3.02652686149083,
      1.405665430359793e-16, 0.2500112482565204,   -0.8942014658868807,  0.1222034745688248,   0.2184605681421391,
      -1.051954866458787,    -0.01451613926015979, -0.6053053722981661,  -3.02652686149083,    1.405665430359793e-16,
      -0.6125634419282462,   -3.06281720964123,    0.1786881422641237,   0.2184605681421391,   -1.051954866458787,
      -0.01451613926015979,  0.2500112482565203,   -0.894201465886881,   -0.1222034745688246,  0.2500112482565204,
      -0.8942014658868807,   0.1222034745688248,   0.2184605681421391,   -1.051954866458787,   -0.01451613926015979,
      -0.01920709001514156,  -0.948045133541213,   0.05910211434006233,  -0.05075777012952284, -1.105798534113119,
      0.04858522096860256,   0.2184605681421391,   -1.051954866458787,   -0.01451613926015979, 0.2500112482565203,
      -0.894201465886881,    -0.1222034745688246,  -0.05075777012952284, -1.105798534113119,   0.04858522096860256,
      -0.01920709001514165,  -0.9480451335412133,  -0.05910211434006224};
}

auto testFile(const std::string& fileName, Dune::Vtk::DataTypes precision, bool native) {
  TestSuite t("Test read in grid");
  auto epsilon = precision == Dune::Vtk::DataTypes::FLOAT32 ? 1e-4 : 1e-8;

  using Grid = Dune::UGGrid<2>;
  Dune::GridFactory<Grid> factory;

  Dune::Vtk::VtkReader reader{factory};
  reader.read(fileName);
  auto grid     = reader.createGrid();
  auto gridView = grid->leafGridView();

  t.check(gridView.size(0) == 4) << testLocation() << "GridView should have 4 elements, but has " << gridView.size(0)
                                 << " elements";

  // get point Data
  auto stressData = reader.getPointData("linearStress");

  t.check(stressData.numComponents() == 3)
      << testLocation() << "Num components should be 3, but is " << stressData.numComponents();
  t.check(stressData.dataType() == precision) << testLocation() << "Precision is wrong";

  auto localStressFunction = localFunction(stressData);
  auto results             = expectedValues(precision, native);
  for (int i = 0; const auto& ele : Dune::elements(gridView)) {
    localStressFunction.bind(ele);
    auto refEle = Dune::referenceElement(ele);
    for (const auto c : Dune::range(ele.geometry().corners())) {
      auto corner       = refEle.geometry<2>(c).center();
      auto stressInFile = localStressFunction(corner);
      for (const auto j : Dune::range(3)) {
        t.check(Dune::FloatCmp::eq(stressInFile[j], results[i], epsilon))
            << "Result read in is " << stressInFile[j] << " but should be " << results[i];
        ++i;
      }
    }
  }

  return t;
}

auto testUnstructuredInstantiaionAndDeduction() {
  TestSuite t("Test ResultFunction");
  std::string fileName = "ResultFunctionTest";

  using Grid     = Dune::UGGrid<2>;
  using GridView = Grid::LeafGridView;

  DummyProblem<Grid> testCase{
      {2, 2}
  };

  auto& gridView        = testCase.gridView();
  auto& sparseAssembler = testCase.sparseAssembler();
  auto& req             = testCase.requirement();
  auto& D_Glob          = req.globalSolution();
  auto fes              = &(testCase.finiteElements());

  // If we initialize result function as VTKFunction we can specify precision, this doesn't affect it when we use it
  // with Dune::Vtk
  auto stressFunction32 =
      Ikarus::makeResultFunction<Ikarus::ResultTypes::linearStress>(sparseAssembler, Dune::VTK::Precision::float32);
  auto stressFunction64 =
      Ikarus::makeResultFunction<Ikarus::ResultTypes::linearStress>(sparseAssembler, Dune::VTK::Precision::float64);
  auto stressVtkFunction = Ikarus::makeResultVtkFunction<Ikarus::ResultTypes::linearStress>(sparseAssembler);

  // Test
  // Dune::VTKWriter will always default to FLOAT32
  Dune::VTKWriter vtkWriter(gridView, Dune::VTK::nonconforming);
  vtkWriter.addVertexData(stressFunction32);
  vtkWriter.write(fileName + "_Common_VTK_32");

  std::cout << "Testing now Common_VTK_32" << std::endl;
  t.subTest(testFile(fileName + "_Common_VTK_32" + ".vtu", Dune::Vtk::DataTypes::FLOAT32, true));

  // You can change coord precision to float64 in the constructor
  Dune::VTKWriter vtkWriter4(gridView, Dune::VTK::nonconforming, Dune::VTK::Precision::float64);
  vtkWriter4.addVertexData(stressFunction64);
  vtkWriter4.write(fileName + "_Common_VTK_64");

  std::cout << "Testing now Common_VTK_64" << std::endl;
  t.subTest(testFile(fileName + "_Common_VTK_64" + ".vtu", Dune::Vtk::DataTypes::FLOAT64, true));

  // Vtk::VtkWriter can be told to use a special Datatype
  Dune::Vtk::UnstructuredGridWriter<GridView, Dune::Vtk::DiscontinuousDataCollector<GridView>> vtkWriter2(
      gridView, Dune::Vtk::FormatTypes::ASCII, Dune::Vtk::DataTypes::FLOAT32);
  vtkWriter2.addPointData(stressVtkFunction);
  auto vtkfileName = vtkWriter2.write(fileName + "_Vtk_32");

  std::cout << "Testing now Vtk_32" << std::endl;
  t.subTest(testFile(vtkfileName, Dune::Vtk::DataTypes::FLOAT32, false));

  Dune::Vtk::UnstructuredGridWriter vtkWriter3(gridView, Dune::Vtk::FormatTypes::ASCII, Dune::Vtk::DataTypes::FLOAT64);
  vtkWriter3.addPointData(stressVtkFunction);
  auto vtkfileName3 = vtkWriter3.write(fileName + "_Vtk_64");

  std::cout << "Testing now Vtk_64" << std::endl;
  t.subTest(testFile(vtkfileName3, Dune::Vtk::DataTypes::FLOAT64, false));

  return t;
}

int main(const int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testUnstructuredInstantiaionAndDeduction());

  return t.exit();
}