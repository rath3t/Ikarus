//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <fstream>
#include <vector>

#include <dune/geometry/type.hh>

#include <Eigen/Core>

#include <ikarus/FiniteElements/ElasticityFE.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Geometries/GeometryType.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>

using namespace Ikarus;

class TestFE {
public:
  static void initialize() {}
  [[nodiscard]] static Ikarus::FiniteElements::IFiniteElement::DofPairVectorType getEntityVariablePairs() {
    return Ikarus::FiniteElements::IFiniteElement::DofPairVectorType{};
  }

  static double calculateScalar(const Ikarus::FiniteElements::ScalarAffordances&, Ikarus::FiniteElements::FEValues&,
                                std::optional<std::reference_wrapper<FiniteElements::FEValues>>& data) {
    if (data)
      if (isType(data->get().get(Ikarus::EntityType::vertex)[0], Ikarus::Variable::VariableTags::displacement2d))
        return getValue(data->get().get(Ikarus::EntityType::vertex)[0])[0]
               * getValue(data->get().get(Ikarus::EntityType::vertex)[0])[1];
    return 5;
  }
};

class TestFE2 {
public:
  static void initialize() {}
  [[nodiscard]] static Ikarus::FiniteElements::IFiniteElement::DofPairVectorType getEntityVariablePairs() {
    return Ikarus::FiniteElements::IFiniteElement::DofPairVectorType{};
  }

  static double calculateScalar(const FiniteElements::ScalarAffordances&, FiniteElements::FEValues&) { return 5; }
};

TEST(FiniteElementInterfaceTest, createGenericFEList) {
  using namespace Ikarus::Grid;
  using namespace Ikarus::FiniteElements;

  using Grid = SimpleGrid<2, 2>;
  SimpleGridFactory<2, 2> gridFactory;
  //  using vertexType = Eigen::Vector2d;
  std::vector<Eigen::Vector2d> verticesVec;
  verticesVec.emplace_back(0.0, 0.0);  // 0
  verticesVec.emplace_back(2.0, 0.0);  // 1
  verticesVec.emplace_back(0.0, 2.0);  // 2
  verticesVec.emplace_back(2.0, 2.0);  // 3
  verticesVec.emplace_back(4.0, 0.0);  // 4
  verticesVec.emplace_back(4.0, 2.0);  // 5

  for (auto&& vert : verticesVec)
    gridFactory.insertVertex(vert);

  std::vector<size_t> elementIndices;
  elementIndices.resize(4);
  elementIndices = {0, 1, 2, 3};
  gridFactory.insertElement(Ikarus::GeometryType::linearQuadrilateral, elementIndices);
  elementIndices = {1, 4, 3, 5};
  gridFactory.insertElement(Ikarus::GeometryType::linearQuadrilateral, elementIndices);

  Grid grid = gridFactory.createGrid();

  auto gridView = grid.leafGridView();
  std::vector<Ikarus::FiniteElements::IFiniteElement> fes;

  for (auto&& element : surfaces(gridView))
    fes.emplace_back(Ikarus::FiniteElements::ElasticityFE(element, gridView.indexSet(), 1000, 0.3));

  Eigen::VectorXd fint{};
  Eigen::MatrixXd K{};
  fint.setZero(8);
  K.setZero(8, 8);
  for (auto&& fe : fes) {
    using namespace Ikarus::Variable;
    std::vector<IVariable> vars;
    vars.emplace_back(VariableFactory::createVariable(VariableTags::displacement2d));
    vars.emplace_back(VariableFactory::createVariable(VariableTags::displacement2d));
    vars.emplace_back(VariableFactory::createVariable(VariableTags::displacement2d));
    vars.emplace_back(VariableFactory::createVariable(VariableTags::displacement2d));

    FiniteElements::FEValues feValues;
    feValues.add(Ikarus::EntityType::vertex, vars);
    feValues.add(Ikarus::EntityType::edge, vars);
    feValues.add(Ikarus::EntityType::surface, vars);
    feValues.add(Ikarus::EntityType::volume, vars);

    // test FE withoutData
    {
      const auto [KEle, fintEle] = calculateLocalSystem(fe, stiffness, forces, feValues);
      EXPECT_EQ(dofSize(fe), 8);
      EXPECT_EQ(calculateVector(fe, forces, feValues).size(), 8);
      EXPECT_DOUBLE_EQ(calculateScalar(fe, potentialEnergy, feValues), 13.0);
      EXPECT_EQ(calculateMatrix(fe, stiffness, feValues).cols(), 8);
      EXPECT_EQ(calculateMatrix(fe, stiffness, feValues).rows(), 8);
      EXPECT_THROW(calculateMatrix(fe, mass, feValues), std::logic_error);
      EXPECT_THROW(calculateLocalSystem(fe, mass, forces, feValues), std::logic_error);
      EXPECT_EQ(KEle.rows(), 8);
      EXPECT_EQ(KEle.cols(), 8);
      EXPECT_EQ(fintEle.size(), 8);
    }
    feValues.add(Ikarus::EntityType::vertex, vars[0]);
    feValues.add(Ikarus::EntityType::edge, vars[0]);
    feValues.add(Ikarus::EntityType::surface, vars[0]);
    feValues.add(Ikarus::EntityType::volume, vars[0]);

    FiniteElements::FEValues dataFeValues;
    dataFeValues.add(Ikarus::EntityType::vertex, vars);
    dataFeValues.add(Ikarus::EntityType::edge, vars);
    dataFeValues.add(Ikarus::EntityType::surface, vars);
    dataFeValues.add(Ikarus::EntityType::volume, vars);

    {
      const auto [KEle, fintEle] = calculateLocalSystem(fe, stiffness, forces, feValues, dataFeValues);
      EXPECT_EQ(dofSize(fe), 8);
      EXPECT_EQ(calculateVector(fe, forces, feValues, dataFeValues).size(), 8);
      EXPECT_DOUBLE_EQ(calculateScalar(fe, potentialEnergy, feValues, dataFeValues), 13.0);
      EXPECT_EQ(calculateMatrix(fe, stiffness, feValues, dataFeValues).cols(), 8);
      EXPECT_EQ(calculateMatrix(fe, stiffness, feValues, dataFeValues).rows(), 8);
      EXPECT_THROW(calculateMatrix(fe, mass, feValues, dataFeValues), std::logic_error);
      EXPECT_THROW(calculateLocalSystem(fe, mass, forces, feValues, dataFeValues), std::logic_error);
      EXPECT_EQ(KEle.rows(), 8);
      EXPECT_EQ(KEle.cols(), 8);
      EXPECT_EQ(fintEle.size(), 8);
    }
  }

  const auto entityIDDofPair = getEntityVariableTuple(fes[0]);
  using namespace Ikarus::Variable;
  std::vector<std::pair<size_t, Ikarus::Variable::VariableTags>> idtagExpected;
  idtagExpected.emplace_back(0, VariableTags::displacement2d);
  idtagExpected.emplace_back(1, VariableTags::displacement2d);
  idtagExpected.emplace_back(2, VariableTags::displacement2d);
  idtagExpected.emplace_back(3, VariableTags::displacement2d);
  for (int i = 0; auto&& [entityID, entityType, varVec] : entityIDDofPair) {
    EXPECT_EQ(entityID, idtagExpected[i].first);
    EXPECT_EQ(varVec.size(), 1);
    EXPECT_EQ(entityType, Ikarus::EntityType::vertex);
    EXPECT_EQ(varVec[0], idtagExpected[i].second);
    ++i;
  }

  auto feT{fes[0]};  // test copy assignment

  using namespace Ikarus::Variable;
  std::vector<IVariable> vars;
  vars.emplace_back(VariableFactory::createVariable(VariableTags::displacement2d));
  vars.emplace_back(VariableFactory::createVariable(VariableTags::displacement2d));
  vars.emplace_back(VariableFactory::createVariable(VariableTags::displacement2d));
  vars.emplace_back(VariableFactory::createVariable(VariableTags::displacement2d));

  FiniteElements::FEValues feValues;
  feValues.add(Ikarus::EntityType::vertex, vars);

  std::vector<IVariable> datas;
  datas.emplace_back(VariableFactory::createVariable(VariableTags::displacement2d));
  datas[0] += Eigen::Vector2d(15, 2);

  FiniteElements::FEValues feDataValues;
  feDataValues.add(Ikarus::EntityType::vertex, datas);

  EXPECT_THAT(feDataValues.get(Ikarus::EntityType::vertex).size(), 1);
  EXPECT_THAT(feDataValues.get(Ikarus::EntityType::edge).size(), 0);
  EXPECT_THAT(feDataValues.get(Ikarus::EntityType::surface).size(), 0);
  EXPECT_THAT(feDataValues.get(Ikarus::EntityType::volume).size(), 0);

  Ikarus::FiniteElements::IFiniteElement fe((TestFE()));
  // check behaviour of dummy fe calculate scalar function before adding data and after
  EXPECT_DOUBLE_EQ(calculateScalar(fe, potentialEnergy, feValues, feDataValues), 30.0);
  datas[0] = VariableFactory::createVariable(VariableTags::displacement1d);
  EXPECT_DOUBLE_EQ(calculateScalar(fe, potentialEnergy, feValues, feDataValues), 5.0);

  Ikarus::FiniteElements::IFiniteElement fe2((TestFE2()));  // check if element without optional data is accepted

  EXPECT_DOUBLE_EQ(calculateScalar(fe2, potentialEnergy, feValues), 5.0);
}