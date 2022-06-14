//
// Created by Alex on 21.07.2021.
//

#include <config.h>

//#include "src/include/ikarus/utils/algorithms.hh"

#include <dune/alugrid/grid.hh>
#include <dune/common/parametertreeparser.hh>
//#include <dune/fufem/dunepython.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <spdlog/spdlog.h>

#include <Eigen/Core>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlRoutines/loadControl.hh>
#include <ikarus/finiteElements/mechanics/3dbeam.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/nonLinearSolver/trustRegion.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/functionSanityChecks.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <ikarus/utils/observer/genericControlObserver.hh>
#include <ikarus/utils/observer/nonLinearSolverLogger.hh>
#include <ikarus/utils/algorithms.hh>

constexpr int centerLineOrder    = 1;
constexpr int quaternionOrder        = 1;
constexpr int gridDim               = 1;
constexpr int worldDim               = 3;
constexpr int quaternionDim           = 4;
constexpr int quaternionCorrectionDim = quaternionDim - 1;

using CenterLinePositionsVector = Dune::BlockVector<Ikarus::RealTuple<double, worldDim>>;
using UnitQuaternionVector      = Dune::BlockVector<Ikarus::UnitVector<double, quaternionDim>>;
using MultiTypeVector = Dune::MultiTypeBlockVector<CenterLinePositionsVector,UnitQuaternionVector >;

int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);
//  Dune::ParameterTree parameterSet;
//  Dune::ParameterTreeParser::readINITree(argv[1], parameterSet);
//
//  const Dune::ParameterTree &gridParameters     = parameterSet.sub("GridParameters");
//  const Dune::ParameterTree &controlParameters  = parameterSet.sub("ControlParameters");
//  const Dune::ParameterTree &materialParameters = parameterSet.sub("MaterialParameters");
//
//  const auto refinement      = gridParameters.get<int>("refinement");
//  const auto innerRadius     = gridParameters.get<double>("innerRadius");
//  const auto mshfilepath     = gridParameters.get<std::string>("mshfilepath");
//  const auto loadSteps       = controlParameters.get<int>("loadSteps");
  const auto loadSteps       = 1;
  const double distrutedLoad       = 10;
//  const auto loadFactorRange = controlParameters.get<std::array<double, 2>>("loadFactorRange");
  const auto loadFactorRange = std::array<double, 2>({0,1});
//
//  const auto A  = materialParameters.get<double>("A");
//  const auto K  = materialParameters.get<double>("K");
//  const auto ms = materialParameters.get<double>("ms");

//  Python::start();
//  Python::Reference main = Python::import("__main__");
//  Python::run("import math");

  using namespace Ikarus;
  Ikarus::BeamMaterial mat({.E = 2.1e8, .nu = 0.0, .A = 0.01, .I1 = 8.333333333*1e-6, .I2 =  8.333333333*1e-6, .J =  16.66666666666*1e-6});

  Dune::GridFactory<Dune::FoamGrid<1, 3, double>> gridFactory;
  const double L = 12.0;
  gridFactory.insertVertex({0, 0,0});
  gridFactory.insertVertex({L,0,0});
  gridFactory.insertElement(Dune::GeometryTypes::line, {0, 1});
  auto grid     = gridFactory.createGrid();
  grid->globalRefine(3);



//  grid->globalRefine(refinement);
  auto gridView = grid->leafGridView();

  //  draw(gridView);
  //  spdlog::info("The domain has a length of {}.", sizedom1);

  using namespace Dune::Functions::BasisFactory;
  //  auto basisEmbedded = makeBasis(gridView, power<directorDim>(lagrange<magnetizationOrder>(),
  //  BlockedInterleaved()));
  auto basisEmbeddedC
      = makeBasis(gridView, composite(power<worldDim>(lagrange<centerLineOrder>(), BlockedInterleaved()),
                                      power<quaternionDim>(lagrange<quaternionOrder>(), BlockedInterleaved()),
                                      BlockedLexicographic{}));

  //  auto basisRie = makeBasis(gridView, power<quaternionCorrectionDim>(lagrange<magnetizationOrder>(),
  //  FlatInterleaved()));
  auto basisRieC = makeBasis(
      gridView, composite(power<worldDim>(lagrange<centerLineOrder>(), FlatInterleaved()),
                          power<quaternionCorrectionDim>(lagrange<quaternionOrder>(), FlatInterleaved()), FlatLexicographic{}));
  std::cout << "This gridview contains: " << std::endl;
  std::cout << gridView.size(1) << " vertices " << std::endl;
  std::cout << gridView.size(0) << " edges" << std::endl;
  std::cout << basisRieC.size() << " Dofs" << std::endl;

  //  draw(gridView);

  std::vector<Ikarus::SimoReissnerBeam<decltype(basisEmbeddedC), decltype(basisRieC)>> fes;
  auto volumeLoad = [distrutedLoad](auto &globalCoord, auto &lamb) {
    Eigen::Vector<double, worldDim> fext;
    fext.setZero();
    fext[0] = 0;
    fext[1] = 0;
    fext[2] = -lamb*distrutedLoad;
    return fext;
  };

  int insideCounter = 0;
  for (auto &element : elements(gridView)) {
    auto geoCoord       = element.geometry().center();
//    std::cout<<geoCoord<<std::endl;
    fes.emplace_back(basisEmbeddedC, basisRieC, element, mat, volumeLoad);
  }

  CenterLinePositionsVector centerLineBlocked(basisEmbeddedC.size({Dune::Indices::_0}));
  for (auto &asingle : centerLineBlocked) {
    asingle.setValue(Eigen::Vector<double, worldDim>::Zero());
  }

  UnitQuaternionVector quatsBlocked(basisEmbeddedC.size({Dune::Indices::_1}));
  for (auto &msingle : quatsBlocked) {
      Eigen::Vector<double, quaternionDim> q;
      q.setZero();
      q.w()=1;
      q.normalize();
      msingle.setValue(q);
  }

  MultiTypeVector mAndABlocked( centerLineBlocked,quatsBlocked);

  std::vector<bool> dirichletFlags(basisRieC.size(), false);
  std::cout << "dirichletFlags.size()" << dirichletFlags.size() << std::endl;
  // Fix displacement on the whole boundary
  Dune::Functions::forEachBoundaryDOF(Dune::Functions::subspaceBasis(basisRieC, Dune::Indices::_0),
                                      [&](auto &&globalIndex) { dirichletFlags[globalIndex[0]] = true; });

    Dune::Functions::forEachBoundaryDOF(Dune::Functions::subspaceBasis(basisRieC, Dune::Indices::_1,1),
                                        [&](auto &&globalIndex) { dirichletFlags[globalIndex[0]] = true; });
   //partly cmaplíng does not work and full clamp also does not work
    Dune::Functions::forEachBoundaryDOF(Dune::Functions::subspaceBasis(basisRieC, Dune::Indices::_1,0),
                                        [&](auto &&globalIndex) { dirichletFlags[globalIndex[0]] = true; });

  // Fix displacement on the whole boundary to zero
  Dune::Functions::forEachBoundaryDOF(
      Dune::Functions::subspaceBasis(basisEmbeddedC, Dune::Indices::_0), [&](auto &&globalIndex) {
        mAndABlocked[Dune::Indices::_0][globalIndex[1]].setValue(Eigen::Vector<double, worldDim>::Zero());
      });

    const auto fixedDofs = std::ranges::count(dirichletFlags, true);

    std::cout << fixedDofs << " FixedDofs" << std::endl;
  auto denseAssembler  = DenseFlatAssembler(basisRieC, fes, dirichletFlags);
  auto sparseAssembler = SparseFlatAssembler(basisRieC, fes, dirichletFlags);
  double lambda        = 0.0;

  auto residualFunction = [&](auto &&disp, auto &&lambdaLocal) -> auto & {
    Ikarus::FErequirements req = FErequirementsBuilder<MultiTypeVector>()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacementAndQuaternions, disp)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::VectorAffordances::forces)
                                     .build();
//      std::cout<<"denseAssembler.getReducedVector(req)"<<std::endl;
//      std::cout<<denseAssembler.getReducedVector(req)<<std::endl;
    return denseAssembler.getReducedVector(req);
  };

  auto hessianFunction = [&](auto &&disp, auto &&lambdaLocal) -> auto & {
    Ikarus::FErequirements req = FErequirementsBuilder<MultiTypeVector>()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacementAndQuaternions, disp)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::MatrixAffordances::stiffness)
                                     .build();
//    std::cout<<"sparseAssembler.getReducedMatrix(req)"<<std::endl;
//    std::cout<<sparseAssembler.getReducedMatrix(req)<<std::endl;
    return sparseAssembler.getReducedMatrix(req);
  };

  auto energyFunction = [&](auto &&disp, auto &&lambdaLocal) -> auto {
    Ikarus::FErequirements req = FErequirementsBuilder<MultiTypeVector>()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacementAndQuaternions, disp)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::ScalarAffordances::mechanicalPotentialEnergy)
                                     .build();

    return denseAssembler.getScalar(req);
  };

  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(energyFunction, residualFunction, hessianFunction),
                                            parameter(mAndABlocked, lambda));
  auto updateFunction = std::function([&](MultiTypeVector &multiTypeVector, const Eigen::VectorXd &d) {
    auto dFull = denseAssembler.createFullVector(d);
    multiTypeVector += dFull;
  });

  checkGradient(nonLinOp, {.draw = false, .writeSlopeStatement = true}, updateFunction);
  checkHessian(nonLinOp, {.draw = false, .writeSlopeStatement = true}, updateFunction);

//  auto nr = Ikarus::makeTrustRegion(nonLinOp, updateFunction);
    auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::sd_UmfPackLU);

    auto nr = Ikarus::makeNewtonRaphson(nonLinOp.subOperator<1, 2>(), std::move(linSolver),updateFunction);

    auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
    nr->subscribeAll(nonLinearSolverObserver);
//    auto nr = Ikarus::makeTrustRegion< decltype(nonLinOp),PreConditioner::DiagonalPreconditioner>(nonLinOp,
//    updateFunction);
//  nr->setup({.verbosity = 1,
//             .maxiter   = 300,
//             .grad_tol  = 1e-7,
//             .corr_tol  = 1e-16,
//             .useRand   = false,
//             .rho_reg   = 1e3,
//             .Delta0    = 1000});

  auto lc = Ikarus::LoadControl(nr, loadSteps, loadFactorRange);

  auto scalarMagnBasis          = makeBasis(gridView, lagrangeDG<centerLineOrder>());
  auto localViewScalarMagnBasis = scalarMagnBasis.localView();

  std::vector<double> gradMNodalRes(scalarMagnBasis.size());
  std::vector<double> divAfieldNodalRes(scalarMagnBasis.size());
  std::vector<Dune::FieldVector<double, 3>> BfieldNodalRes(scalarMagnBasis.size());
  std::vector<Dune::FieldVector<double, 3>> HfieldNodalRes(scalarMagnBasis.size());

  auto writerObserver = std::make_shared<Ikarus::GenericControlObserver>(ControlMessages::STEP_ENDED, [&](auto i) {
    auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, worldDim>>(
        Dune::Functions::subspaceBasis(basisEmbeddedC, Dune::Indices::_0), mAndABlocked);
    auto quatFunc = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, quaternionDim>>(
        Dune::Functions::subspaceBasis(basisEmbeddedC, Dune::Indices::_1), mAndABlocked);

    Dune::VTKWriter vtkWriter(gridView, Dune::VTK::conforming);
    vtkWriter.addVertexData(disp, Dune::VTK::FieldInfo("m", Dune::VTK::FieldInfo::Type::vector, worldDim));
//    vtkWriter.addVertexData(quatFunc, Dune::VTK::FieldInfo("A", Dune::VTK::FieldInfo::Type::vector, vectorPotDim));

//    auto resultRequirements
//        = Ikarus::ResultRequirementsBuilder<MultiTypeVector>()
//              .insertGlobalSolution(Ikarus::FESolutions::displacementAndQuaternions, mAndABlocked)
//              .insertParameter(Ikarus::FEParameter::loadfactor, lambda)
//              .addResultRequest(ResultType::gradientNormOfMagnetization, ResultType::BField, ResultType::HField,
//                                ResultType::divergenceOfVectorPotential)
//              .build();
//    auto localmFunction = localFunction(mGlobalFunc);
//
//    auto ele = elements(gridView).begin();
//    ResultTypeMap<double> result;
////    for (auto &fe : fes) {
////      localViewScalarMagnBasis.bind(*ele);
////      const auto &fe2              = localViewScalarMagnBasis.tree().finiteElement();
////      const auto &referenceElement = Dune::ReferenceElements<double, gridDim>::general(ele->type());
////      for (auto c = 0UL; c < fe2.size(); ++c) {
////        const auto fineKey                        = fe2.localCoefficients().localKey(c);
////        const auto nodalPositionInChildCoordinate = referenceElement.position(fineKey.subEntity(), fineKey.codim());
////
////        auto coord = toEigenVector(nodalPositionInChildCoordinate);
////        fe.calculateAt(resultRequirements, coord, result);
////        gradMNodalRes[localViewScalarMagnBasis.index(localViewScalarMagnBasis.tree().localIndex(c))[0]]
////            = result.getResult(ResultType::gradientNormOfMagnetization)(0, 0);
////        BfieldNodalRes[localViewScalarMagnBasis.index(localViewScalarMagnBasis.tree().localIndex(c))[0]]
////            = toFieldVector(result.getResult(ResultType::BField));
////        HfieldNodalRes[localViewScalarMagnBasis.index(localViewScalarMagnBasis.tree().localIndex(c))[0]]
////            = toFieldVector(result.getResult(ResultType::HField));
////        divAfieldNodalRes[localViewScalarMagnBasis.index(localViewScalarMagnBasis.tree().localIndex(c))[0]]
////            = result.getResult(ResultType::divergenceOfVectorPotential)(0, 0);
////      }
////      ++ele;
////    }
//
//    auto gradmGlobalFunc = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(scalarMagnBasis, gradMNodalRes);
//    auto divAGlobalFunc  = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(scalarMagnBasis, divAfieldNodalRes);
//    auto bFieldGlobalFunc = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(
//        scalarMagnBasis, BfieldNodalRes);
//    auto hFieldGlobalFunc = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(
//        scalarMagnBasis, HfieldNodalRes);
//
//    vtkWriter.addVertexData(gradmGlobalFunc, Dune::VTK::FieldInfo("gradMNorm", Dune::VTK::FieldInfo::Type::scalar, 1));
//    vtkWriter.addVertexData(bFieldGlobalFunc, Dune::VTK::FieldInfo("B", Dune::VTK::FieldInfo::Type::vector, 3));
//    vtkWriter.addVertexData(hFieldGlobalFunc, Dune::VTK::FieldInfo("H", Dune::VTK::FieldInfo::Type::vector, 3));
//    vtkWriter.addVertexData(divAGlobalFunc, Dune::VTK::FieldInfo("divA", Dune::VTK::FieldInfo::Type::scalar, 1));
//    auto isInsideFunc = Dune::Functions::makeAnalyticGridViewFunction(isInsidePredicate, gridView);
//    vtkWriter.addCellData(isInsideFunc, Dune::VTK::FieldInfo("isInside", Dune::VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.write(std::string("Beam") + std::to_string(i));

    std::cout << "Writing " << std::string("Beam") + std::to_string(i) << " to file. The load factor is " << lambda
              << std::endl;
  });
  lc.subscribeAll(writerObserver);
  lc.run();
  nonLinOp.update<0>();
  spdlog::info("Energy: {}", nonLinOp.value());
  spdlog::info("MaxDisplacement: {}", std::ranges::max(nonLinOp.firstParameter()[Dune::Indices::_0],[](auto& e1,auto& e2) {return std::abs(e1[2])<std::abs(e2[2]);})[2]);
    std::cout<<mAndABlocked<<std::endl;
    checkGradient(nonLinOp, {.draw = false, .writeSlopeStatement = true}, updateFunction);
    checkHessian(nonLinOp, {.draw = false, .writeSlopeStatement = true}, updateFunction);
    spdlog::info("Max Disp Bernoulli beam: {}", 5.0/384.0*distrutedLoad*L*L*L*L/(mat.E*mat.I1));
}