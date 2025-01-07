// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <ikarus/utils/functionhelper.hh>

using Dune::TestSuite;

#include <dune/foamgrid/foamgrid.hh>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mechanics/loads/volume.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>

template <typename TestSuiteType, typename SparseType, typename DenseType>
void checkAssembledQuantities(TestSuiteType& t, SparseType& sType, DenseType& dType, std::size_t dofSize) {
  t.check(isApproxSame(sType, dType, 1e-15), "Dense==Sparse");
  t.check(sType.rows() == static_cast<Eigen::Index>(dofSize))
      << "DOFsCheck via rows: " << sType.rows() << "rows and " << dofSize << " DOFs";
  if (not(std::is_same_v<std::remove_cvref_t<SparseType>, Eigen::VectorXd>))
    t.check(sType.cols() == static_cast<Eigen::Index>(dofSize))
        << "DOFsCheck via columns: " << sType.cols() << " cols and " << dofSize << " DOFs";
}

template <typename PreBasis>
auto SimpleAssemblersTest(const PreBasis& preBasis) {
  TestSuite t("SimpleAssemblersTest");
  using namespace Ikarus;
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox       = {4, 2};
  std::array<int, 2> elementsPerDirection = {2, 1};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);

  for (int ref = 0; ref < 4; ++ref) {
    auto gridView = grid->leafGridView();

    auto basis     = Ikarus::makeBasis(gridView, preBasis);
    auto totalDOFs = basis.flat().size();

    auto vL = []([[maybe_unused]] const auto& globalCoord, const auto& lamb) {
      Eigen::Vector2d fext;
      fext.setZero();
      fext[1] = 2 * lamb;
      fext[0] = lamb;
      return fext;
    };

    auto linMat = Ikarus::LinearElasticity(Ikarus::toLamesFirstParameterAndShearModulus({.emodul = 100, .nu = 0.2}));
    auto sk     = skills(linearElastic(planeStress(linMat)), volumeLoad<2>(vL));

    using FEType = decltype(makeFE(basis, sk));
    std::cout << Dune::className<FEType>() << std::endl;

    std::vector<FEType> fes;
    for (auto&& ge : elements(gridView)) {
      fes.emplace_back(makeFE(basis, sk));
      fes.back().bind(ge);
    }

    using ChildType = std::remove_cvref_t<decltype(fes[0].localView().tree().child(0))>;

    auto basisP = std::make_shared<const decltype(basis)>(basis);
    Ikarus::DirichletValues dirichletValues(basisP->flat());
    dirichletValues.fixDOFs([](auto& basis_, auto& dirichletFlags) {
      Dune::Functions::forEachBoundaryDOF(basis_, [&](auto&& indexGlobal) { dirichletFlags[indexGlobal] = true; });
    });

    // center position to be fixed in all directions
    Dune::FieldVector<double, 2> pos{2.0, 1.0};
    int centerNode = 1; // number of nodes at the center ( = pos)
    if (Ikarus::Concepts::LagrangeNodeOfOrder<ChildType, 1> and ref == 0) {
      t.checkThrow<Dune::InvalidStateException>(
          [&]() { auto fixIndices = utils::globalIndexFromGlobalPosition(basis.flat(), pos); },
          "globalIndexFromGlobalPosition should have failed for order = 1 and ref = 0 as no node exists at "
          "the center.");
      centerNode = 0;
    } else {
      const auto fixIndices = utils::globalIndexFromGlobalPosition(basis.flat(), pos);
      for (const auto idx : fixIndices)
        dirichletValues.setSingleDOF(idx, true);
    }

    using SparseAssembler = SparseFlatAssembler<decltype(fes), decltype(dirichletValues)>;
    using ScalarAssembler = ScalarFlatAssembler<decltype(fes), decltype(dirichletValues)>;
    using VectorAssembler = VectorFlatAssembler<decltype(fes), decltype(dirichletValues)>;
    using DenseAssembler  = DenseFlatAssembler<decltype(fes), decltype(dirichletValues)>;

    SparseAssembler sparseFlatAssembler(fes, dirichletValues);
    ScalarAssembler scalarFlatAssembler(fes, dirichletValues);
    VectorAssembler vectorFlatAssembler(fes, dirichletValues);
    DenseAssembler denseFlatAssembler(fes, dirichletValues);

    auto sparseFlatAssemblerAM = makeAssemblerManipulator(sparseFlatAssembler);
    auto scalarFlatAssemblerAM = makeAssemblerManipulator(scalarFlatAssembler);
    auto vectorFlatAssemblerAM = makeAssemblerManipulator(vectorFlatAssembler);
    auto denseFlatAssemblerAM  = makeAssemblerManipulator(denseFlatAssembler);

    // AssemblerManipulator<ScalarAssembler> scalarFlatAssemblerAM(fes, dirichletValues);
    // AssemblerManipulator<VectorAssembler> vectorFlatAssemblerAM(fes, dirichletValues);
    // AssemblerManipulator<DenseAssembler> denseFlatAssemblerAM(fes, dirichletValues);

    Eigen::VectorXd d(basis.flat().size());
    d.setRandom();
    double load = 0.0;

    auto req = typename FEType::Requirement(d, load);

    const auto& KRawDense = denseFlatAssembler.matrix(req, MatrixAffordance::stiffness, DBCOption::Raw);
    const auto& KRaw      = sparseFlatAssembler.matrix(req, MatrixAffordance::stiffness, DBCOption::Raw);
    const auto& RRawDense = denseFlatAssembler.vector(req, VectorAffordance::forces, DBCOption::Raw);
    const auto& RRaw      = sparseFlatAssembler.vector(req, VectorAffordance::forces, DBCOption::Raw);
    const auto& RVecRaw   = vectorFlatAssembler.vector(req, VectorAffordance::forces, DBCOption::Raw);
    checkAssembledQuantities(t, KRaw, KRawDense, totalDOFs);
    checkAssembledQuantities(t, RRaw, RRawDense, totalDOFs);
    checkAssembledQuantities(t, RVecRaw, RRawDense, totalDOFs);

    const auto& KDense = denseFlatAssembler.matrix(req, MatrixAffordance::stiffness, DBCOption::Full);
    const auto& K      = sparseFlatAssembler.matrix(req, MatrixAffordance::stiffness, DBCOption::Full);
    const auto& RDense = denseFlatAssembler.vector(req, VectorAffordance::forces, DBCOption::Full);
    const auto& R      = sparseFlatAssembler.vector(req, VectorAffordance::forces, DBCOption::Full);
    const auto& RVec   = vectorFlatAssembler.vector(req, VectorAffordance::forces, DBCOption::Full);
    checkAssembledQuantities(t, K, KDense, totalDOFs);
    checkAssembledQuantities(t, R, RDense, totalDOFs);
    checkAssembledQuantities(t, RVec, RDense, totalDOFs);
    const double energyDense  = denseFlatAssembler.scalar(req, ScalarAffordance::mechanicalPotentialEnergy);
    const double energy       = sparseFlatAssembler.scalar(req, ScalarAffordance::mechanicalPotentialEnergy);
    const double energyScalar = scalarFlatAssembler.scalar(req, ScalarAffordance::mechanicalPotentialEnergy);
    checkScalars(t, energy, energyDense, " Incorrect energies for sparse and dense assemblers");
    checkScalars(t, energyScalar, energyDense, " Incorrect energies for scalar and dense assemblers");

    const Eigen::Index fixedDOFs = dirichletValues.fixedDOFsize();
    int boundaryNodes            = (elementsPerDirection[0] * Dune::power(2, ref) + 1) * 2 +
                        (elementsPerDirection[1] * Dune::power(2, ref) + 1) * 2 - 4;
    if constexpr (Ikarus::Concepts::LagrangeNodeOfOrder<ChildType, 2>)
      boundaryNodes *= 2;

    t.check(2 * (boundaryNodes + centerNode) == fixedDOFs) << "Boundary DOFs (" << 2 * (boundaryNodes + centerNode)
                                                           << ") is not equal to Fixed DOFs (" << fixedDOFs << ")";

    /// check if full matrices and full vectors are correct after applying boundary conditions
    t.check(std::ranges::count(KDense.reshaped(), 1) == fixedDOFs) << "Correct number of ones in matrix";
    t.check(std::ranges::count(RDense, 0) == fixedDOFs) << "Correct number of zeros in vector";

    for (auto i = 0U; i < totalDOFs; ++i) {
      if (dirichletValues.isConstrained(i)) {
        for (auto j = 0U; j < totalDOFs; ++j) {
          t.check(Dune::FloatCmp::eq(KDense(i, j), static_cast<double>(i == j)))
              << std::string(i == j ? "" : "Off-") + "Diagonal term is " << KDense(i, j) << " at (" << i << "," << j
              << ")";
        }
        t.check(Dune::FloatCmp::eq(RDense(i), 0.0)) << "Vector component is " << RDense(i) << " at " << i;
      } else {
        for (auto j = 0U; j < totalDOFs; ++j)
          if (not(dirichletValues.isConstrained(j)))
            t.check(Dune::FloatCmp::eq(KDense(i, j), KRawDense(i, j)))
                << "KDense and KRawDense have the same entries corresponding to unconstrained DOFs";
        t.check(Dune::FloatCmp::eq(RDense(i), RRawDense(i)))
            << "RDense and RRawDense have the same entries corresponding to unconstrained DOFs";
      }
    }

    const auto& KRedDense = denseFlatAssembler.matrix(req, MatrixAffordance::stiffness, DBCOption::Reduced);
    const auto& KRed      = sparseFlatAssembler.matrix(req, MatrixAffordance::stiffness, DBCOption::Reduced);
    const auto& RRedDense = denseFlatAssembler.vector(req, VectorAffordance::forces, DBCOption::Reduced);
    const auto& RRed      = sparseFlatAssembler.vector(req, VectorAffordance::forces, DBCOption::Reduced);
    const auto& RVecRed   = vectorFlatAssembler.vector(req, VectorAffordance::forces, DBCOption::Reduced);
    checkAssembledQuantities(t, KRed, KRedDense, totalDOFs - fixedDOFs);
    checkAssembledQuantities(t, RRed, RRedDense, totalDOFs - fixedDOFs);
    checkAssembledQuantities(t, RVecRed, RRedDense, totalDOFs - fixedDOFs);

    /// check if reduced matrices and reduced vectors are correct after applying boundary conditions
    Eigen::Index r = 0U;
    for (auto i = 0U; i < totalDOFs; ++i) {
      if (not(dirichletValues.isConstrained(i))) {
        Eigen::Index c = 0U;
        for (auto j = 0U; j < totalDOFs; ++j)
          if (not(dirichletValues.isConstrained(j))) {
            t.check(Dune::FloatCmp::eq(KRawDense(i, j), KRedDense(r, c)))
                << "Matrix components are " << KRawDense(i, j) << " and " << KRedDense(r, c) << " at (" << r << "," << c
                << ")";
            c++;
          }
        t.check(Dune::FloatCmp::eq(RRawDense(i), RRedDense(r)))
            << "Vector components are " << RRawDense(i) << " and " << RRedDense(r) << " at " << r;
        r++;
      }
    }

    if (not(Ikarus::Concepts::LagrangeNodeOfOrder<ChildType, 1> and ref == 0)) {
      const auto globalIndices = utils::globalIndexFromGlobalPosition(basis.flat(), pos);
      double springStiffness   = 596482;
      double pointF            = 64.23;
      auto doubleEnergy        = [&](const auto&, const auto&, auto, double& energyL) -> void { energyL *= 2.0; };
      auto pointLoad           = [&](const auto&, const auto&, auto, auto, Eigen::VectorXd& vec) -> void {
        vec[globalIndices[0]] += pointF;
        vec[globalIndices[1]] += pointF;
      };
      auto addSpringStiffnessDense = [&](const auto&, const auto&, auto, auto, Eigen::MatrixXd& mat) -> void {
        mat.diagonal()[globalIndices[0]] += springStiffness;
        mat.diagonal()[globalIndices[1]] += springStiffness;
      };
      auto addSpringStiffnessSparse = [&](const auto&, const auto&, auto, auto,
                                          Eigen::SparseMatrix<double>& mat) -> void {
        mat.diagonal()[globalIndices[0]] += springStiffness;
        mat.diagonal()[globalIndices[1]] += springStiffness;
      };

      scalarFlatAssemblerAM->bind(doubleEnergy);

      vectorFlatAssemblerAM->bind(pointLoad);

      sparseFlatAssemblerAM->bind(doubleEnergy);
      sparseFlatAssemblerAM->bind(pointLoad);
      sparseFlatAssemblerAM->bind(addSpringStiffnessSparse);

      denseFlatAssemblerAM->bind(doubleEnergy);
      denseFlatAssemblerAM->bind(doubleEnergy); // doubled twice
      denseFlatAssemblerAM->bind(pointLoad);
      denseFlatAssemblerAM->bind(addSpringStiffnessDense);

      const auto& mKRawDense = denseFlatAssemblerAM->matrix(req, MatrixAffordance::stiffness, DBCOption::Raw);
      const auto& mKRaw      = sparseFlatAssemblerAM->matrix(req, MatrixAffordance::stiffness, DBCOption::Raw);
      const auto& mRRawDense = denseFlatAssemblerAM->vector(req, VectorAffordance::forces, DBCOption::Raw);
      const auto& mRRaw      = sparseFlatAssemblerAM->vector(req, VectorAffordance::forces, DBCOption::Raw);
      const auto& mRVecRaw   = vectorFlatAssemblerAM->vector(req, VectorAffordance::forces, DBCOption::Raw);
      checkAssembledQuantities(t, mKRaw, mKRawDense, totalDOFs);
      checkAssembledQuantities(t, mRRaw, mRRawDense, totalDOFs);
      checkAssembledQuantities(t, mRVecRaw, mRRawDense, totalDOFs);

      const auto& mKDense = denseFlatAssemblerAM->matrix(req, MatrixAffordance::stiffness, DBCOption::Full);
      const auto& mK      = sparseFlatAssemblerAM->matrix(req, MatrixAffordance::stiffness, DBCOption::Full);
      const auto& mRDense = denseFlatAssemblerAM->vector(req, VectorAffordance::forces, DBCOption::Full);
      const auto& mR      = sparseFlatAssemblerAM->vector(req, VectorAffordance::forces, DBCOption::Full);
      const auto& mRVec   = vectorFlatAssemblerAM->vector(req, VectorAffordance::forces, DBCOption::Full);
      checkAssembledQuantities(t, mK, mKDense, totalDOFs);
      checkAssembledQuantities(t, mR, mRDense, totalDOFs);
      checkAssembledQuantities(t, mRVec, mRDense, totalDOFs);

      const double mEnergyDense = denseFlatAssemblerAM->scalar(req, ScalarAffordance::mechanicalPotentialEnergy);
      const double mEnergy      = sparseFlatAssemblerAM->scalar(req, ScalarAffordance::mechanicalPotentialEnergy);
      const double mEnergySca   = scalarFlatAssemblerAM->scalar(req, ScalarAffordance::mechanicalPotentialEnergy);
      checkScalars(t, mEnergy * 2, mEnergyDense, " Incorrect energies for manipulated sparse and dense assemblers");
      checkScalars(t, mEnergySca * 2, mEnergyDense, " Incorrect energies for manipulated scalar and dense assemblers");
      checkScalars(t, mEnergy, 2.0 * energy, " Incorrect energy for sparse assembler after manipulation");
      checkScalars(t, mEnergyDense, 4.0 * energyDense, " Incorrect energy for dense assembler after manipulation");
      for (int i = 0; i < 2; ++i) {
        auto index = globalIndices[i];
        checkScalars(t, KRawDense.diagonal()[index] + springStiffness, mKRawDense.diagonal()[index],
                     " Mismatch in KRawDense and mKRawDense at index = " + std::to_string(i));
        checkScalars(t, KRaw.diagonal()[index] + springStiffness, mKRaw.diagonal()[index],
                     " Mismatch in KRaw and mKRaw at index = " + std::to_string(i));
        checkScalars(t, KDense.diagonal()[index] + springStiffness, mKDense.diagonal()[index],
                     " Mismatch in KDense and mKDense at index = " + std::to_string(i));
        checkScalars(t, K.diagonal()[index] + springStiffness, mK.diagonal()[index],
                     " Mismatch in K and mK at index = " + std::to_string(i));
        checkScalars(t, RRawDense[index] + pointF, mRRawDense[index],
                     " Mismatch in RRawDense and mRRawDense at index = " + std::to_string(i));
        checkScalars(t, RRaw[index] + pointF, mRRaw[index],
                     " Mismatch in RRaw and mRRaw at index = " + std::to_string(i));
        checkScalars(t, RVecRaw[index] + pointF, mRVecRaw[index],
                     " Mismatch in RVecRaw and mRVecRaw at index = " + std::to_string(i));
        checkScalars(t, RDense[index] + pointF, mRDense[index],
                     " Mismatch in RDense and mRDense at index = " + std::to_string(i));
        checkScalars(t, R[index] + pointF, mR[index], " Mismatch in R and mR at index = " + std::to_string(i));
        checkScalars(t, RVec[index] + pointF, mRVec[index],
                     " Mismatch in RVec and mRVec at index = " + std::to_string(i));
      }
    }

    grid->globalRefine(1);
  }
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  using namespace Dune::Functions::BasisFactory;
  auto firstOrderLagrangePrePower2Basis  = power<2>(lagrange<1>(), FlatInterleaved());
  auto secondOrderLagrangePrePower2Basis = power<2>(lagrange<2>(), FlatInterleaved());

  t.subTest(SimpleAssemblersTest(firstOrderLagrangePrePower2Basis));
  t.subTest(SimpleAssemblersTest(secondOrderLagrangePrePower2Basis));
  return t.exit();
}
