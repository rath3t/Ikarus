// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <vector>

#include <dune/alugrid/grid.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/common/typetraits.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/vtk/vtkwriter.hh>

#include "ikarus/assembler/simpleassemblers.hh"
#include "ikarus/utils/dirichletvalues.hh"

#if HAVE_DUNE_IGA
  #include <dune/iga/nurbsgrid.hh>
#endif
#if HAVE_DUNE_LOCALFEFUNCTIONS
  #include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
#endif
#include "testhelpers.hh"

#include <ikarus/finiteelements/autodifffe.hh>
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
struct OneDFoamGridIn2D
{
};
struct OneDFoamGridIn3D
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
  } else if constexpr (std::is_same_v<GridType, Grids::OneDFoamGridIn2D>) {
    Dune::GridFactory<Dune::FoamGrid<1, 2, double>> gridFactory;
    constexpr double h = 1.0;
    constexpr double L = 10.0;
    gridFactory.insertVertex({0, 0});
    gridFactory.insertVertex({L, h});
    gridFactory.insertVertex({2 * L, 0});
    gridFactory.insertElement(Dune::GeometryTypes::line, {0, 1});
    gridFactory.insertElement(Dune::GeometryTypes::line, {1, 2});
    auto grid = gridFactory.createGrid();
    return grid;
  } else if constexpr (std::is_same_v<GridType, Grids::OneDFoamGridIn3D>) {
    Dune::GridFactory<Dune::FoamGrid<1, 3, double>> gridFactory;
    gridFactory.insertVertex({0, 0, 0});
    gridFactory.insertVertex({2, 1, 4});
    gridFactory.insertVertex({4, 0, 1});
    gridFactory.insertElement(Dune::GeometryTypes::line, {0, 1});
    gridFactory.insertElement(Dune::GeometryTypes::line, {1, 2});
    auto grid = gridFactory.createGrid();
    return grid;
  } else
    DUNE_THROW(Dune::NotImplemented, "The requested GridType is not implemented.");
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
      << "calculateVector is not the gradient of calculateScalar." << messageIfFailed;
  return t;
}

template <typename NonLinearOperator>
[[nodiscard]] auto checkHessianOfElement(NonLinearOperator& nonLinearOperator,
                                         const std::string& messageIfFailed = "") {
  Dune::TestSuite t("Check Hessian");
  t.check(Ikarus::utils::checkHessian(nonLinearOperator, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << "calculateMatrix is not the Hessian of calculateScalar. " << messageIfFailed;
  return t;
}

template <typename NonLinearOperator>
[[nodiscard]] auto checkJacobianOfElement(NonLinearOperator& nonLinearOperator,
                                          const std::string& messageIfFailed = "") {
  Dune::TestSuite t("Check Jacobian");
  t.check(Ikarus::utils::checkJacobian(nonLinearOperator, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << "The Jacobian of calculateVector is not calculateMatrix." << messageIfFailed;
  return t;
}

template <template <typename, int, int> class RT, bool vectorizedResult = true>
[[nodiscard]] auto checkCalculateAt(auto& /*nonLinearOperator*/, auto& fe, const auto& feRequirements,
                                    const auto& expectedResult, const auto& evaluationPositions,
                                    const std::string& messageIfFailed = "") {
  Dune::TestSuite t("Test of the calulateAt function for " + Dune::className(fe), Dune::TestSuite::AlwaysThrow);

  using FiniteElement = std::remove_cvref_t<decltype(fe)>;
  Eigen::MatrixXd computedResults(expectedResult.rows(), expectedResult.cols());

  if constexpr (requires { fe.template calculateAt<RT>(feRequirements, evaluationPositions[0]); }) {
    for (int i = 0; const auto& pos : evaluationPositions) {
      if constexpr (vectorizedResult) {
        auto result              = fe.template calculateAt<RT>(feRequirements, pos).asVec();
        computedResults.row(i++) = result.transpose();
      } else {
        auto result              = fe.template calculateAt<RT>(feRequirements, pos).asMat();
        computedResults.row(i++) = result.reshaped().transpose();
      }
    }
    const bool isResultCorrect = isApproxSame(computedResults, expectedResult, 1e-8);
    t.check(isResultCorrect) << "Computed Result for " << Ikarus::toString<RT>()
                             << " is not the same as expected result:\n"
                             << "It is:\n"
                             << computedResults << "\nBut should be:\n"
                             << expectedResult << "\n"
                             << messageIfFailed;
  } else
    static_assert(Dune::AlwaysFalse<FiniteElement>::value, "Element can not provide the requested ResultType ");

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

  // Make a dummy assembler that holds the elements
  Ikarus::DirichletValues dirichletValues(fe.localView().globalBasis());
  auto sparseAssembler = Ikarus::makeSparseFlatAssembler(fes, dirichletValues);
  sparseAssembler->bind(feRequirements);
  sparseAssembler->bind(Ikarus::DBCOption::Full);

  auto vtkResultFunction   = Ikarus::makeResultVtkFunction<resType, ResultEvaluator>(sparseAssembler);
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
  auto resultFunction = Ikarus::makeResultFunction<resType, ResultEvaluator>(sparseAssembler);

  vtkWriter2.addVertexData(resultFunction);
  vtkWriter2.write("native_vtkwriter_resultfunction_" + resultFunction->name() + "_" +
                   std::to_string(FiniteElement::myDim) + std::to_string(element.geometry().type().id()));

  return t;
}

template <typename NonLinearOperator, typename FiniteElement,
          typename FERequirementType = typename FiniteElement::FERequirementType, typename AffordanceColl>
[[nodiscard]] auto checkFEByAutoDiff(NonLinearOperator&, FiniteElement& fe, FERequirementType req,
                                     AffordanceColl affordance, const std::string& messageIfFailed = "") {
  Dune::TestSuite t("Check calculateScalarImpl() and calculateVectorImpl() by Automatic Differentiation");
  auto& basis           = fe.localView().globalBasis();
  auto nDOF             = basis.size();
  using AutoDiffBasedFE = Ikarus::AutoDiffFE<FiniteElement, true>;
  AutoDiffBasedFE feAutoDiff(fe);

  const double tol = 1e-10;

  Eigen::MatrixXd K, KAutoDiff;
  K.setZero(nDOF, nDOF);
  KAutoDiff.setZero(nDOF, nDOF);

  calculateMatrix(fe, req, affordance.matrixAffordance(), K);
  calculateMatrix(feAutoDiff, req, affordance.matrixAffordance(), KAutoDiff);

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

    calculateVector(fe, req, affordance.vectorAffordance(), R);
    calculateVector(feAutoDiff, req, affordance.vectorAffordance(), RAutoDiff);
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
    auto energy         = calculateScalar(fe, req, affordance.scalarAffordance());
    auto energyAutoDiff = calculateScalar(feAutoDiff, req, affordance.scalarAffordance());
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
