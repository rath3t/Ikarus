// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <Eigen/Core>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/functionhelper.hh>

template <typename Material>
constexpr bool isPlaneStress =
    (Ikarus::traits::isSpecializationNonTypeAndTypes<Ikarus::VanishingStress, Material>::value);

inline auto linearStressResultsOfSquare = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 4;
  constexpr int quantities = 3;

  Eigen::Matrix<double, vertices, quantities> expectedStress;
  if constexpr (requires { fe.setEASType(4); }) {
    fe.setEASType(4);
    if (isPlaneStress<typename FE::Material>)
      expectedStress << 1214.28571429, 1214.28571429, 384.61538462, 1214.28571429, 214.28571429, 384.61538462,
          214.28571429, 1214.28571429, 384.61538462, 214.28571429, 214.28571429, 384.61538462;
    else
      expectedStress << 1510.98901099, 1510.98901099, 384.61538462, 1510.98901099, 412.08791209, 384.61538462,
          412.08791209, 1510.98901099, 384.61538462, 412.08791209, 412.08791209, 384.61538462;
  } else {
    if (isPlaneStress<typename FE::Material>)
      expectedStress << 1428.57142857, 1428.57142857, 769.23076923, 1098.90109890, 329.67032967, 384.61538462,
          329.67032967, 1098.90109890, 384.61538462, 0, 0, 0;
    else
      expectedStress << 1923.07692308, 1923.07692308, 769.23076923, 1346.15384615, 576.92307692, 384.61538462,
          576.92307692, 1346.15384615, 384.61538462, 0, 0, 0;
  }

  auto& displacement = nonLinearOperator.firstParameter();
  displacement << 0, 0, 1, 1, 1, 1, 1, 1;

  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearVonMisesResultsOfSquare = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 4;
  constexpr int quantities = 1;

  const auto expectedStress =
      isPlaneStress<typename FE::Material>
          ? Eigen::Matrix<double, vertices, quantities>{1953.44932249, 1182.27663689, 1182.27663689, 0}
          : Eigen::Matrix<double, vertices, quantities>{2339.52405011, 1346.15384615, 1346.15384615, 0};

  auto& displacement = nonLinearOperator.firstParameter();
  displacement << 0, 0, 1, 1, 1, 1, 1, 1;

  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearPrincipalStressResultsOfSquare = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 4;
  constexpr int quantities = 2;

  const auto expectedStress = isPlaneStress<typename FE::Material>
   ? Eigen::Matrix<double, vertices, quantities>{
      {2197.80219780, 659.34065934},
      {1258.21400751, 170.35742107},
      {1258.21400751, 170.35742107},
      {            0,            0}
  } 
  : Eigen::Matrix<double, vertices, quantities>{
    {2692.30769231, 1153.84615385},
    {1505.46675476, 417.61016832},
    {1505.46675476, 417.61016832},
    {0, 0}
  };

  auto& displacement = nonLinearOperator.firstParameter();
  displacement << 0, 0, 1, 1, 1, 1, 1, 1;

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

  auto& displacement = nonLinearOperator.firstParameter();
  displacement << 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearVonMisesResultsOfCube = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 8;
  constexpr int quantities = 1;

  const Eigen::Matrix<double, vertices, quantities> expectedStress{1216.26063853, 0, 2035.19331620, 1216.26063853, 0, 0,
                                                                   1216.26063853, 0};

  auto& displacement = nonLinearOperator.firstParameter();
  displacement << 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

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

  auto& displacement = nonLinearOperator.firstParameter();
  displacement << 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearStressResultsOfTriangle = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 3;
  constexpr int quantities = 3;

  Eigen::Matrix<double, vertices, quantities> expectedStress;
  if (isPlaneStress<typename FE::Material>)
    expectedStress.rowwise() = Eigen::Matrix<double, 1, quantities>{2197.80219780, 659.34065934, 384.61538462};
  else
    expectedStress.rowwise() = Eigen::Matrix<double, 1, quantities>{2692.30769231, 1153.84615385, 384.61538462};

  auto& displacement = nonLinearOperator.firstParameter();
  displacement << 0, 0, 2, 0, 1, 0;

  auto feRequirements = typename FE::Requirement().insertGlobalSolution(displacement);

  return std::make_tuple(feRequirements, expectedStress, Ikarus::utils::referenceElementVertexPositions(fe));
};

inline auto linearStressResultsOfTetrahedron = []<typename NOP, typename FE>(NOP& nonLinearOperator, FE& fe) {
  constexpr int vertices   = 4;
  constexpr int quantities = 6;

  Eigen::Matrix<double, vertices, quantities> expectedStress;
  expectedStress.rowwise() =
      Eigen::Matrix<double, 1, quantities>{576.92307692, 1346.15384615, 576.92307692, 0, 384.61538462, 769.23076923};

  auto& displacement = nonLinearOperator.firstParameter();
  displacement << 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0;

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