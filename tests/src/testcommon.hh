// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <vector>

#include <dune/alugrid/grid.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/vtk/vtkwriter.hh>
#if HAVE_DUNE_IGA
  #include <dune/iga/nurbsgrid.hh>
#endif
#if HAVE_DUNE_LOCALFEFUNCTIONS
  #include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
#endif
#include "testhelpers.hh"

#include <ikarus/finiteelements/autodiff/autodifffe.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/feresulttypes.hh>
#include <ikarus/io/resultfunction.hh>
#include <ikarus/utils/eigendunetransformations.hh>
#include <ikarus/utils/functionsanitychecks.hh>
#include <ikarus/utils/linearalgebrahelper.hh>
#include <ikarus/utils/pythonautodiffdefinitions.hh>

namespace Grids {
struct Yasp
{
};
struct Alu
{
};
struct IgaSurfaceIn2D
{
};
struct IgaSurfaceIn3D
{
};
} // namespace Grids

template <typename GridType>
auto createGrid([[maybe_unused]] int elex = 10, [[maybe_unused]] int eley = 10) {
  // //  /// ALUGrid Example
  if constexpr (std::is_same_v<GridType, Grids::Alu>) {
    using Grid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
    auto grid  = Dune::GmshReader<Grid>::read("testfiles/unstructuredtrianglesfine.msh", false);
    grid->globalRefine(0);
    return grid;
  } else if constexpr (std::is_same_v<GridType, Grids::Yasp>) {
    using Grid     = Dune::YaspGrid<2>;
    const double L = 1;
    const double h = 1;

    Dune::FieldVector<double, 2> bbox       = {L, h};
    std::array<int, 2> elementsPerDirection = {elex, eley};
    auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);
    return grid;
  } else if constexpr (std::is_same_v<GridType, Grids::IgaSurfaceIn2D>) {
#if HAVE_DUNE_IGA
    constexpr auto dimworld = 2;
    const std::array order  = {2, 2};

    const std::array<std::vector<double>, 2> knotSpans = {
        {{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}
    };

    using ControlPoint = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointType;

    const std::vector<std::vector<ControlPoint>> controlPoints = {
        {  {.p = {0, 0}, .w = 5},    {.p = {0.5, 0}, .w = 1},   {.p = {1, 0}, .w = 1}},
        {{.p = {0, 0.5}, .w = 1}, {.p = {0.5, 0.5}, .w = 10}, {.p = {1, 0.5}, .w = 1}},
        {  {.p = {0, 1}, .w = 1},    {.p = {0.5, 1}, .w = 1},   {.p = {1, 1}, .w = 1}}
    };

    std::array<int, 2> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};

    auto controlNet = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointNetType(dimsize, controlPoints);
    using Grid      = Dune::IGA::NURBSGrid<2, dimworld>;

    Dune::IGA::NURBSPatchData<2, dimworld> patchData;
    patchData.knotSpans     = knotSpans;
    patchData.degree        = order;
    patchData.controlPoints = controlNet;
    auto grid               = std::make_shared<Grid>(patchData);
    grid->globalRefine(1);
    return grid;
  } else if constexpr (std::is_same_v<GridType, Grids::IgaSurfaceIn3D>) {
    constexpr auto dimworld = 3;
    const std::array order  = {2, 2};

    const std::array<std::vector<double>, 2> knotSpans = {
        {{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}
    };

    using ControlPoint = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointType;

    const std::vector<std::vector<ControlPoint>> controlPoints = {
        {   {.p = {0, 0, 0.1}, .w = 5},   {.p = {0.5, 0, 0.1}, .w = 1},   {.p = {1, 0, 0.2}, .w = 1}},
        {{.p = {0, 0.5, -0.2}, .w = 1},  {.p = {0.5, 0.5, 1}, .w = 10}, {.p = {1, 0.5, 0.2}, .w = 1}},
        {   {.p = {0, 1, 0.4}, .w = 1}, {.p = {0.5, 1, -0.25}, .w = 1},   {.p = {1, 1, 0.7}, .w = 1}}
    };

    std::array<int, 2> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};

    auto controlNet = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointNetType(dimsize, controlPoints);
    using Grid      = Dune::IGA::NURBSGrid<2, dimworld>;

    Dune::IGA::NURBSPatchData<2, dimworld> patchData;
    patchData.knotSpans     = knotSpans;
    patchData.degree        = order;
    patchData.controlPoints = controlNet;
    auto grid               = std::make_shared<Grid>(patchData);
    grid->globalRefine(1);
    return grid;
#endif
  }
}

template <int size>
struct CornerFactory
{
  static void construct(std::vector<Dune::FieldVector<double, size>>& values, const int corners = 10) {
    values.resize(corners);
    std::generate(values.begin(), values.end(),
                  []() { return Ikarus::createRandomVector<Dune::FieldVector<double, size>>(); });
  }
};

enum class CornerDistortionFlag
{
  unDistorted,
  fixedDistorted,
  randomlyDistorted
};

// Corner factory for element with codim==0, e.g. no surfaces in 3D
template <int gridDim>
struct ValidCornerFactory
{
  static void construct(std::vector<Dune::FieldVector<double, gridDim>>& values, const Dune::GeometryType& type,
                        const CornerDistortionFlag& distortionFlag) {
    const auto& refElement = Dune::ReferenceElements<double, gridDim>::general(type);

    const auto numberOfVertices = refElement.size(gridDim);

    values.resize(numberOfVertices);
    for (int i = 0; i < numberOfVertices; ++i)
      values[i] = refElement.position(i, gridDim);

    /// Perturbation of the corner values
    if (distortionFlag == CornerDistortionFlag::unDistorted)
      return;
    else if (distortionFlag == CornerDistortionFlag::fixedDistorted) {
      if (not((gridDim == 2) and (numberOfVertices == 4)))
        DUNE_THROW(Dune::NotImplemented, "Fixed distortion is only implemented for a 2D 4-node element (Q1)");
      std::vector<Dune::FieldVector<double, gridDim>> randomnessOnNodes;
      randomnessOnNodes.push_back({-0.2, -0.05});
      randomnessOnNodes.push_back({-0.15, 0.05});
      randomnessOnNodes.push_back({0.15, 0.15});
      randomnessOnNodes.push_back({-0.05, -0.1});
      for (long i = 0; i < 4; ++i)
        values[i] += randomnessOnNodes[i];
    } else if (distortionFlag == CornerDistortionFlag::randomlyDistorted) {
      std::transform(values.begin(), values.end(), values.begin(), [](const auto& vec) {
        return vec + Ikarus::createRandomVector<Dune::FieldVector<double, gridDim>>(-0.2, 0.2);
      });
    } else
      DUNE_THROW(Dune::IOError, "Incorrect CornerDistortionFlag entered");
  }
};

template <int gridDim>
auto createUGGridFromCorners(const CornerDistortionFlag& distortionFlag,
                             const Dune::GeometryType& geometryType = Dune::GeometryTypes::cube(gridDim)) {
  using Grid = Dune::UGGrid<gridDim>;

  std::vector<Dune::FieldVector<double, gridDim>> corners;

  const int numberOfVertices = Dune::ReferenceElements<double, gridDim>::general(geometryType).size(gridDim);

  ValidCornerFactory<gridDim>::construct(corners, geometryType, distortionFlag);

  std::vector<unsigned int> vertexArrangment;
  vertexArrangment.resize(numberOfVertices);
  std::iota(vertexArrangment.begin(), vertexArrangment.end(), 0);

  Dune::GridFactory<Grid> gridFactory;
  for (auto& corner : corners) {
    gridFactory.insertVertex(corner);
  }
  gridFactory.insertElement(geometryType, vertexArrangment);

  std::unique_ptr<Grid> grid = gridFactory.createGrid();
  return grid;
}

template <typename Ele>
struct ElementTest
{
};

template <typename NonLinearOperator>
[[nodiscard]] auto checkGradientOfElement(NonLinearOperator& nonLinearOperator,
                                          const std::string& messageIfFailed = "") {
  Dune::TestSuite t("Check gradient");
  t.check(Ikarus::utils::checkGradient(nonLinearOperator, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << "The gradient of calculateVector is not the gradient of calculateScalar." << messageIfFailed;
  return t;
}

template <typename NonLinearOperator>
[[nodiscard]] auto checkHessianOfElement(NonLinearOperator& nonLinearOperator,
                                         const std::string& messageIfFailed = "") {
  Dune::TestSuite t("Check Hessian");
  t.check(Ikarus::utils::checkHessian(nonLinearOperator, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << "The Hessian of calculateMatrix is not the Hessian of calculateScalar. " << messageIfFailed;
  return t;
}

template <typename NonLinearOperator>
[[nodiscard]] auto checkJacobianOfElement(NonLinearOperator& nonLinearOperator,
                                          const std::string& messageIfFailed = "") {
  Dune::TestSuite t("Check Jacobian");
  t.check(Ikarus::utils::checkJacobian(nonLinearOperator, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << "The Jacobian of calculateMatrix is not the Jacobian of calculateVector." << messageIfFailed;
  return t;
}

template <template <typename, int, int> class resType, bool voigt = true>
[[nodiscard]] auto checkCalculateAt(auto& /*nonLinearOperator*/, auto& fe, const auto& feRequirements,
                                    const auto& expectedResult, const auto& evaluationPositions,
                                    const std::string& messageIfFailed = "") {
  Dune::TestSuite t("Test of the calulateAt function for " + Dune::className(fe));

  using FiniteElement = std::remove_cvref_t<decltype(fe)>;
  Eigen::MatrixXd computedResults(expectedResult.rows(), expectedResult.cols());

  if constexpr (FiniteElement::template canProvideResultType<resType>()) {
    for (int i = 0; const auto& pos : evaluationPositions) {
      if constexpr (voigt) {
        auto result              = fe.template calculateAt<resType>(feRequirements, pos).asVec();
        computedResults.row(i++) = result.transpose();
      } else {
        auto result = fe.template calculateAt<resType>(feRequirements, pos).asMat();
        Eigen::Map<Eigen::VectorXd> vector(result.data(), result.size());
        computedResults.row(i++) = vector.transpose();
      }
    }
    const bool isResultCorrect = isApproxSame(computedResults, expectedResult, 1e-8);
    t.check(isResultCorrect) << "Computed Result for " << Ikarus::toString<resType>()
                             << " is not the same as expected result:\n"
                             << "It is:\n"
                             << computedResults << "\nBut should be:\n"
                             << expectedResult << "\n"
                             << messageIfFailed;
  } else
    t.check(false) << "Element can not provide the requested RsultType " << Ikarus::toString<resType>()
                   << messageIfFailed;

  return t;
}

template <template <typename, int, int> class resType, typename ResultEvaluator>
[[nodiscard]] auto checkResultFunction(auto& /*nonLinearOperator*/, auto& fe, const auto& feRequirements,
                                       const auto& expectedResult, const auto& evaluationPositions,
                                       const std::string& messageIfFailed = "") {
  Dune::TestSuite t("Result Function Test" + Dune::className(fe));

  using FiniteElement = std::remove_reference_t<decltype(fe)>;
  auto& element       = fe.gridElement();
  auto gridView       = fe.localView().globalBasis().gridView();
  std::vector<FiniteElement> fes{fe};

  Eigen::MatrixXd computedResults(expectedResult.rows(), expectedResult.cols());

  auto vtkResultFunction   = Ikarus::makeResultVtkFunction<resType, ResultEvaluator>(&fes, feRequirements);
  auto localResultFunction = localFunction(vtkResultFunction);
  localResultFunction.bind(element);

  for (int i = 0; const auto& pos : evaluationPositions) {
    auto result = localResultFunction(pos);
    for (auto j : std::views::iota(0ul, result.size()))
      computedResults(i, j) = result[j];
    ++i;
  }
  const bool isResultCorrect = isApproxSame(computedResults, expectedResult, 1e-8);
  t.check(isResultCorrect) << "Computed Result for " << Ikarus::toString<resType>()
                           << " is not the same as expected result:\n"
                           << "It is:\n"
                           << computedResults << "\nBut should be:\n"
                           << expectedResult << "\n"
                           << messageIfFailed;

  Dune::Vtk::VtkWriter vtkWriter(gridView);

  vtkWriter.addPointData(vtkResultFunction);
  vtkWriter.write("Vtk_VtkWriter_resultfunction_" + vtkResultFunction.name() + "_" +
                  std::to_string(FiniteElement::myDim) + std::to_string(element.geometry().type().id()));

  Dune::VTKWriter<decltype(gridView)> vtkWriter2(gridView);
  auto resultFunction = Ikarus::makeResultFunction<resType, ResultEvaluator>(&fes, feRequirements);

  vtkWriter2.addVertexData(resultFunction);
  vtkWriter2.write("native_vtkwriter_resultfunction_" + resultFunction->name() + "_" +
                   std::to_string(FiniteElement::myDim) + std::to_string(element.geometry().type().id()));

  return t;
}

template <typename NonLinearOperator, typename FiniteElement,
          typename FERequirementType = typename FiniteElement::FERequirementType>
[[nodiscard]] auto checkFEByAutoDiff(NonLinearOperator&, FiniteElement& fe, FERequirementType req,
                                     const std::string& messageIfFailed = "") {
  Dune::TestSuite t("Check calculateScalarImpl() and calculateVectorImpl() by Automatic Differentiation");
  auto& basis           = fe.localView().globalBasis();
  auto nDOF             = basis.size();
  using AutoDiffBasedFE = Ikarus::AutoDiffFE<FiniteElement, true>;
  AutoDiffBasedFE feAutoDiff(fe);

  const double tol = 1e-10;

  Eigen::MatrixXd K, KAutoDiff;
  K.setZero(nDOF, nDOF);
  KAutoDiff.setZero(nDOF, nDOF);

  fe.calculateMatrix(req, K);
  feAutoDiff.calculateMatrix(req, KAutoDiff);

  if constexpr (requires { feAutoDiff.getFE().getNumberOfEASParameters(); }) {
    t.check(fe.getNumberOfEASParameters() == feAutoDiff.getFE().getNumberOfEASParameters())
        << "Number of EAS parameters for FE(" << fe.getNumberOfEASParameters()
        << ") and number of EAS parameters for AutodiffFE(" << feAutoDiff.getFE().getNumberOfEASParameters()
        << ") are not equal";
  }

  t.check(K.isApprox(KAutoDiff, tol),
          "Mismatch between the stiffness matrices obtained from explicit implementation and the one based on "
          "automatic differentiation")
      << messageIfFailed << "\nKAutoDiff:\n"
      << KAutoDiff << "\nK:\n"
      << K;

  try {
    Eigen::VectorXd R, RAutoDiff;
    R.setZero(nDOF);
    RAutoDiff.setZero(nDOF);

    fe.calculateVector(req, R);
    feAutoDiff.calculateVector(req, RAutoDiff);
    t.check(R.isApprox(RAutoDiff, tol),
            "Mismatch between the residual vectors obtained from explicit implementation and the one based on "
            "automatic differentiation")
        << messageIfFailed << "\nRAutoDiff:\n"
        << RAutoDiff << "\nR:\n"
        << R;
  } catch (const Dune::NotImplemented&) {
    /// calculateVector not checked since the element throws a Dune::NotImplemented
  }

  try {
    auto energy         = fe.calculateScalar(req);
    auto energyAutoDiff = feAutoDiff.calculateScalar(req);
    t.check(Dune::FloatCmp::eq(energy, energyAutoDiff, tol),
            "Mismatch between the energies obtained from explicit implementation and the one based on "
            "automatic differentiation")
        << messageIfFailed;
  } catch (const Dune::NotImplemented&) {
  }

  return t;
}

template <typename NonLinearOperator>
void resetNonLinearOperatorParametersToZero(NonLinearOperator& nonLinOp) {
  nonLinOp.firstParameter().setZero();
  nonLinOp.lastParameter() = 0.0;
  nonLinOp.updateAll();
}
