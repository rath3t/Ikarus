// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <Eigen/Core>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/mechanics/materials/vanishingstress.hh>
#include <ikarus/utils/functionhelper.hh>

namespace Testing {

constexpr double NaN = std::numeric_limits<double>::signaling_NaN();
template <typename Material>
constexpr bool isPlaneStress =
    (Ikarus::traits::isSpecializationNonTypeAndTypes<Ikarus::Materials::VanishingStress, Material>::value);

static Eigen::Vector<double, 8> displacementsForSquare({0, 0, 1, 1, 1, 1, 1, 1});
static Eigen::Vector<double, 24> displacementsForCube({0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0,
                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
static Eigen::Vector<double, 6> displacementsForTriangle({0, 0, 2, 0, 1, 0});
static Eigen::Vector<double, 12> displacementsForTetrahedron({0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0});

} // namespace Testing

inline auto linearStressResultsOfSquare = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 4;
  constexpr int quantities = 3;

  Eigen::Matrix<double, vertices, quantities> expectedStress;
  if constexpr (requires { fe.setEASType(4); }) {
    fe.setEASType(4);
    if (Testing::isPlaneStress<typename FE::Material>)
      expectedStress << 1214.28571429, 1214.28571429, 384.61538462, 1214.28571429, 214.28571429, 384.61538462,
          214.28571429, 1214.28571429, 384.61538462, 214.28571429, 214.28571429, 384.61538462;
    else
      expectedStress << 1510.98901099, 1510.98901099, 384.61538462, 1510.98901099, 412.08791209, 384.61538462,
          412.08791209, 1510.98901099, 384.61538462, 412.08791209, 412.08791209, 384.61538462;
  } else {
    if (Testing::isPlaneStress<typename FE::Material>)
      expectedStress << 1428.57142857, 1428.57142857, 769.23076923, 1098.90109890, 329.67032967, 384.61538462,
          329.67032967, 1098.90109890, 384.61538462, 0, 0, 0;
    else
      expectedStress << 1923.07692308, 1923.07692308, 769.23076923, 1346.15384615, 576.92307692, 384.61538462,
          576.92307692, 1346.15384615, 384.61538462, 0, 0, 0;
  }

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForSquare;
  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linear3dPlaneStrainStressResultsOfSquare = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  static_assert(!Testing::isPlaneStress<typename FE::Material>);

  constexpr int vertices   = 4;
  constexpr int quantities = 6;

  Eigen::Matrix<double, vertices, quantities> expectedStress{
      {1923.07692308, 1923.07692308, 1153.84615385, 0, 0, 769.23076923},
      {1346.15384615,  576.92307692,  576.92307692, 0, 0, 384.61538462},
      { 576.92307692, 1346.15384615,  576.92307692, 0, 0, 384.61538462},
      {            0,             0,             0, 0, 0,            0}
  };

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForSquare;
  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearPolarStressResultsOfSquare = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  static_assert(Testing::isPlaneStress<typename FE::Material>);

  constexpr int vertices   = 4;
  constexpr int quantities = 3;

  Eigen::Matrix<double, vertices, quantities> expectedStress{
      {2197.80219780,  659.34065934,   -0.00000000},
      { 329.67032967, 1098.90109890,  384.61538462},
      { 329.67032967, 1098.90109890, -384.61538462},
      {            0,             0,             0},
  };

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForSquare;
  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearVonMisesResultsOfSquare = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 4;
  constexpr int quantities = 1;

  const auto expectedStress =
      Testing::isPlaneStress<typename FE::Material>
          ? Eigen::Matrix<double, vertices, quantities>{1953.44932249, 1182.27663689, 1182.27663689, 0}
          : Eigen::Matrix<double, vertices, quantities>{1538.46153846, 1017.59665810, 1017.59665810, 0};

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForSquare;
  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearHydrostaticStressResultsOfSquare = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 4;
  constexpr int quantities = 1;

  const auto expectedStress =
      Testing::isPlaneStress<typename FE::Material>
          ? Eigen::Matrix<double, vertices, quantities>{1428.57142857, 714.28571429, 714.28571429, 0}
          : Eigen::Matrix<double, vertices, quantities>{1666.66666667, 833.33333333, 833.33333333, 0};

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForSquare;
  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearTriaxialityStressResultsOfSquare = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 4;
  constexpr int quantities = 1;

  const auto expectedStress =
      Testing::isPlaneStress<typename FE::Material>
          ? Eigen::Matrix<double, vertices, quantities>{0.73130714, 0.60416124, 0.60416124, Testing::NaN}
          : Eigen::Matrix<double, vertices, quantities>{1.08333333, 0.81892302, 0.81892302, Testing::NaN};

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForSquare;
  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearPrincipalStressResultsOfSquare = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices = 4;

  const auto expectedStress = []() {
    if constexpr (Testing::isPlaneStress<typename FE::Material>)
      return Eigen::Matrix<double, vertices, 2>{
          {2197.80219780, 659.34065934},
          {1258.21400751, 170.35742107},
          {1258.21400751, 170.35742107},
          {            0,            0}
      };
    else
      return Eigen::Matrix<double, vertices, 3>{
          {2692.30769231, 1153.84615385, 1153.84615385},
          {1505.46675476,  576.92307692,  417.61016832},
          {1505.46675476,  576.92307692,  417.61016832},
          {            0,             0,             0}
      };
  }();

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForSquare;
  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearStressResultsOfCube = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 8;
  constexpr int quantities = 6;

  const Eigen::Matrix<double, vertices, quantities> expectedStress{
      {  576.92307692, 1346.15384615,   576.92307692,  384.61538462,             0,  384.61538462},
      {             0,             0,              0,             0,             0,             0},
      {-1346.15384615,  192.30769231, -1346.15384615,             0, -769.23076923,             0},
      {-1346.15384615, -576.92307692,  -576.92307692,             0, -384.61538462, -384.61538462},
      {             0,             0,              0,             0,             0,             0},
      {             0,             0,              0,             0,             0,             0},
      { -576.92307692, -576.92307692, -1346.15384615, -384.61538462, -384.61538462,             0},
      {             0,             0,              0,             0,             0,             0}
  };

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForCube;
  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearVonMisesResultsOfCube = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 8;
  constexpr int quantities = 1;

  const Eigen::Matrix<double, vertices, quantities> expectedStress{1216.26063853, 0, 2035.19331620, 1216.26063853, 0, 0,
                                                                   1216.26063853, 0};

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForCube;
  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearTriaxialityResultsOfCube = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 8;
  constexpr int quantities = 1;

  const Eigen::Matrix<double, vertices, quantities> expectedStress{
      0.68516016, Testing::NaN, -0.40946151, -0.68516016, Testing::NaN, Testing::NaN, -0.68516016, Testing::NaN};

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForCube;

  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearHydrostaticStressResultsOfCube = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 8;
  constexpr int quantities = 1;

  const Eigen::Matrix<double, vertices, quantities> expectedStress{833.33333333,  0, -833.33333333, -833.33333333, 0, 0,
                                                                   -833.33333333, 0};

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForCube;
  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearPrincipalStressResultsOfCube = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 8;
  constexpr int quantities = 3;

  Eigen::Matrix<double, vertices, quantities> expectedStress{
      {1627.71184906,  576.92307692,   295.36507401},
      {            0,             0,              0},
      { 192.30769231, -576.92307692, -2115.38461538},
      {-295.36507401, -576.92307692, -1627.71184906},
      {            0,             0,              0},
      {            0,             0,              0},
      {-295.36507401, -576.92307692, -1627.71184906},
      {            0,             0,              0}
  };

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForCube;
  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearStressResultsOfTriangle = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 3;
  constexpr int quantities = 3;

  Eigen::Matrix<double, vertices, quantities> expectedStress;
  if (Testing::isPlaneStress<typename FE::Material>)
    expectedStress.rowwise() = Eigen::Matrix<double, 1, quantities>{2197.80219780, 659.34065934, 384.61538462};
  else
    expectedStress.rowwise() = Eigen::Matrix<double, 1, quantities>{2692.30769231, 1153.84615385, 384.61538462};

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForTriangle;

  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearPolarStressResultsOfTriangle = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  static_assert(Testing::isPlaneStress<typename FE::Material>);

  constexpr int vertices   = 3;
  constexpr int quantities = 3;

  Eigen::Matrix<double, vertices, quantities> expectedStress{
      {1813.18681319, 1043.95604396, -769.23076923},
      {1582.41758242, 1274.72527473,  846.15384615},
      { 659.34065934, 2197.80219780,  384.61538462}
  };

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForTriangle;
  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearStressResultsOfTetrahedron = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 4;
  constexpr int quantities = 6;

  Eigen::Matrix<double, vertices, quantities> expectedStress;
  expectedStress.rowwise() =
      Eigen::Matrix<double, 1, quantities>{576.92307692, 1346.15384615, 576.92307692, 0, 384.61538462, 769.23076923};

  typename NOP::Domain::SolutionVectorType displacement = Testing::displacementsForTetrahedron;
  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

template <typename CompileTimeMatrix>
auto stressResultsToMatrix(const CompileTimeMatrix& expectedResults) {
  constexpr int vertices   = CompileTimeMatrix::CompileTimeTraits::RowsAtCompileTime;
  constexpr int voigtComps = CompileTimeMatrix::CompileTimeTraits::ColsAtCompileTime;
  constexpr int matrixSize = (-1 + Ikarus::ct_sqrt(1 + 8 * vertices)) / 2;
  constexpr int comps      = matrixSize * matrixSize;

  Eigen::Matrix<double, vertices, comps> transformedResults;
  for (const auto i : std::views::iota(0, vertices)) {
    auto stressMatrix         = Ikarus::fromVoigt(Eigen::Vector<double, voigtComps>(expectedResults.row(i)), false);
    transformedResults.row(i) = Eigen::Map<Eigen::VectorXd>(stressMatrix.data(), comps);
  }
  return transformedResults;
}