// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/finiteelements/feresulttypes.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/makeenum.hh>

using Dune::TestSuite;
using namespace Ikarus;

MAKE_ENUM(ElementHasResultType, Symmetric, Unsymmetric, Vector, Scalar, Custom, Crazy, Non_Square)

// Matrices
REGISTER_SIMPLE_SYMMETRIC_RESULTTYPE(Res1, 3, 3, false);
REGISTER_SIMPLE_RESULTTYPE(Res2, 3, 3);

// Vectors
REGISTER_SIMPLE_RESULTTYPE(Res3, 3, 1);

// Scalars
REGISTER_SIMPLE_RESULTTYPE(Res4, 1, 1);

// Dynamic
REGISTER_SIMPLE_RESULTTYPE(Res5, Eigen::Dynamic, Eigen::Dynamic);

// Fixed Size Max
REGISTER_RESERVED_RESULTTYPE(Res6, Eigen::Dynamic, Eigen::Dynamic, 7, 13);

// Non-Square matrix
REGISTER_SIMPLE_RESULTTYPE(Res7, 3, 2);

template <ElementHasResultType rt, Ikarus::ResultShape shape>
struct DummyElement
{
  constexpr static bool asVec = shape == Ikarus::ResultShape::Vector;

  static auto calculateAt()
  requires(rt == ElementHasResultType::Symmetric)
  {
    ResultWrapper<Res1<double, 1, 1>, shape> result;
    if constexpr (asVec)
      result = Eigen::Matrix<double, 6, 1>::Zero();
    else
      result = Eigen::Matrix<double, 3, 3>::Zero();
    return result;
  }
  static auto calculateAt()
  requires(rt == ElementHasResultType::Unsymmetric)
  {
    ResultWrapper<Res2<double, 1, 1>, shape> result;
    if constexpr (asVec)
      result = Eigen::Matrix<double, 9, 1>::Zero();
    else
      result = Eigen::Matrix<double, 3, 3>::Zero();
    return result;
  }
  static auto calculateAt()
  requires(rt == ElementHasResultType::Vector)
  {
    ResultWrapper<Res3<double, 1, 1>, shape> result;
    result = Eigen::Matrix<double, 3, 1>::Zero();
    return result;
  }
  static auto calculateAt()
  requires(rt == ElementHasResultType::Scalar)
  {
    ResultWrapper<Res4<double, 1, 1>, shape> result;
    result = Eigen::Matrix<double, 1, 1>::Zero();
    return result;
  }
  static auto calculateAt()
  requires(rt == ElementHasResultType::Custom)
  {
    ResultWrapper<Res5<double, 1, 1>, shape> result;
    if constexpr (asVec)
      result = Eigen::MatrixXd::Zero(9, 1);
    else
      result = Eigen::MatrixXd::Zero(3, 3);
    return result;
  }

  static auto calculateAt()
  requires(rt == ElementHasResultType::Crazy)
  {
    ResultWrapper<Res6<double, 1, 1>, shape> result;
    if constexpr (asVec)
      result = Eigen::MatrixXd::Zero(66, 1);
    else
      result = Eigen::MatrixXd::Zero(6, 11);
    return result;
  }

  static auto calculateAt()
  requires(rt == ElementHasResultType::Non_Square)
  {
    ResultWrapper<Res7<double, 1, 1>, shape> result;
    if constexpr (asVec)
      result = Eigen::Matrix<double, 6, 1>::Zero();
    else
      result = Eigen::Matrix<double, 3, 2>::Zero();
    return result;
  }
};

struct Shape
{
  int rows;
  int cols;
};

auto testRTs() {
  static_assert(Concepts::ResultType<Res1>);
  static_assert(Concepts::ResultType<Res2>);
  static_assert(Concepts::ResultType<Res3>);
  static_assert(Concepts::ResultType<Res4>);
  static_assert(Concepts::ResultType<Res5>);
  static_assert(Concepts::ResultType<Res6>);
  static_assert(Concepts::ResultType<Res7>);
  TestSuite t("Test FE ResultTypes");

  auto testShapes = [&]<ElementHasResultType rt, Ikarus::ResultShape shape>(Shape expectedShapeVec,
                                                                            Shape expectedShapeMat) {
    auto resultContainer = DummyElement<rt, shape>::calculateAt();
    auto resultAsVec     = resultContainer.asVec();
    auto resultAsMat     = resultContainer.asMat(expectedShapeMat.rows, expectedShapeMat.cols);

    TestSuite tSub("Test Case: " + toString(rt));
    tSub.check(expectedShapeVec.rows == resultAsVec.rows())
        << "Expected rows for Vec is " << expectedShapeVec.rows << " but is " << resultAsVec.rows();
    tSub.check(expectedShapeVec.cols == resultAsVec.cols())
        << "Expected cols for Vec is " << expectedShapeVec.cols << " but is " << resultAsVec.cols();

    tSub.check(expectedShapeMat.rows == resultAsMat.rows())
        << "Expected rows for Mat is " << expectedShapeMat.rows << " but is " << resultAsMat.rows();
    tSub.check(expectedShapeMat.cols == resultAsMat.cols())
        << "Expected cols for Mat is " << expectedShapeMat.cols << " but is " << resultAsMat.cols();
    t.subTest(tSub);
  };

  // Symmetric case
  testShapes.operator()<ElementHasResultType::Symmetric, ResultShape::Vector>({6, 1}, {3, 3});
  testShapes.operator()<ElementHasResultType::Symmetric, ResultShape::Matrix>({6, 1}, {3, 3});

  // Unsymmetric case
  testShapes.operator()<ElementHasResultType::Unsymmetric, ResultShape::Vector>({9, 1}, {3, 3});
  testShapes.operator()<ElementHasResultType::Unsymmetric, ResultShape::Matrix>({9, 1}, {3, 3});

  // Vector case
  testShapes.operator()<ElementHasResultType::Vector, ResultShape::Vector>({3, 1}, {3, 1});

  // Scalar case
  testShapes.operator()<ElementHasResultType::Scalar, ResultShape::Vector>({1, 1}, {1, 1});

  // Dynamic case
  testShapes.operator()<ElementHasResultType::Custom, ResultShape::Matrix>({9, 1}, {3, 3});
  testShapes.operator()<ElementHasResultType::Custom, ResultShape::Vector>({9, 1}, {3, 3});

  // Crazy case
  testShapes.operator()<ElementHasResultType::Crazy, ResultShape::Matrix>({66, 1}, {6, 11});
  testShapes.operator()<ElementHasResultType::Crazy, ResultShape::Vector>({66, 1}, {6, 11});

  // Non_Square case
  testShapes.operator()<ElementHasResultType::Non_Square, ResultShape::Matrix>({6, 1}, {3, 2});
  testShapes.operator()<ElementHasResultType::Non_Square, ResultShape::Vector>({6, 1}, {3, 2});

  return t;
}

int main(const int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testRTs());

  return t.exit();
}
