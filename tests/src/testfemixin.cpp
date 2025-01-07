// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>

using Dune::TestSuite;

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/materials/svk.hh>
#include <ikarus/finiteelements/mechanics/nonlinearelastic.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>

template <typename PreBasis>
auto mixinTest(const PreBasis& preBasis) {
  TestSuite t("SimpleAssemblersTest");
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox       = {4, 2};
  std::array<int, 2> elementsPerDirection = {2, 1};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);

  auto gridView = grid->leafGridView();

  auto basis        = Ikarus::makeBasis(gridView, preBasis);
  auto matParameter = Ikarus::toLamesFirstParameterAndShearModulus({.emodul = 1000, .nu = 0.3});
  auto totalDOFs    = basis.flat().size();

  Ikarus::Materials::StVenantKirchhoff matSVK(matParameter);
  auto reducedMat = Ikarus::Materials::planeStress(matSVK, 1e-8);
  using namespace Ikarus;

  auto vL = []([[maybe_unused]] const Dune::FieldVector<double, 2>& globalCoord, const double& lamb) {
    Eigen::Vector2d fext;
    fext.setZero();
    fext[1] = 2 * lamb;
    fext[0] = lamb;
    return fext;
  };

  Dune::BitSetVector<1> neumannVertices(gridView.size(2), true);
  BoundaryPatch neumannBoundary(gridView, neumannVertices);

  auto nBLoad = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, const auto& lamb) {
    Eigen::Vector2d fext;
    fext.setZero();
    fext[0] = lamb / 40;
    return fext;
  };

  auto preFE = makeFE(
      basis, skills(nonLinearElastic(reducedMat), volumeLoad(vL), neumannBoundaryLoad(&neumannBoundary, nBLoad)));

  using FEType = decltype(preFE);
  std::vector<FEType> fes;

  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(preFE);
    fes.back().bind(ge);
  }

  Eigen::VectorXd d(basis.flat().size());

  d.setRandom();
  double load = 0.0;

  auto req = typename FEType::Requirement(d, load);

  Eigen::VectorXd forces(fes.begin()->size());
  forces.setZero();
  Eigen::MatrixXd stiffness(fes.begin()->size(), fes.begin()->size());
  stiffness.setZero();
  const FE firstFE = *fes.begin();
  calculateVector(firstFE, req, VectorAffordance::forces, forces);
  calculateMatrix(firstFE, req, MatrixAffordance::stiffness, stiffness);
  Dune::FieldVector<double, 2> local;
  static_assert(not requires { firstFE.template calculateAt<ResultTypes::linearStress>(req, local); });
  static_assert(requires { firstFE.template calculateAt<ResultTypes::PK2Stress>(req, local); });

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  using namespace Dune::Functions::BasisFactory;
  auto firstOrderLagrangePrePower2Basis  = power<2>(lagrange<1>(), FlatInterleaved());
  auto secondOrderLagrangePrePower2Basis = power<2>(lagrange<2>(), FlatInterleaved());

  t.subTest(mixinTest(firstOrderLagrangePrePower2Basis));
  t.subTest(mixinTest(secondOrderLagrangePrePower2Basis));
  return t.exit();
}
