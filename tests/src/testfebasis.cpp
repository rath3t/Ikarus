// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/nedelecbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/raviartthomasbasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <ikarus/finiteelements/fefactory.hh>

using Dune::TestSuite;

#include <ikarus/finiteelements/febase.hh>
#include <ikarus/finiteelements/fehelper.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>

/**
 * \brief A helper function to get a FEBase.
 *
 * \tparam useFlat A boolean indicating if the type of the underlying basis is of the flat or the untouched version.
 * \tparam GridView Type of the grid view
 * \tparam PreBasis Type of the pre basis
 * *
 * \param gridView The DUNE grid view.
 * \param preBasis The pre basis to make a basis.
 *
 * \return FEBase corresponding to the gridView and preBasis.
 */
template <bool useFlat, typename GridView, typename PreBasis>
auto getFEBase(const GridView& gridView, const PreBasis& preBasis) {
  const auto basis = Ikarus::makeBasis(gridView, preBasis);
  auto element     = elements(gridView).begin();
  auto fe          = Ikarus::makeFE<useFlat>(basis, Ikarus::skills());
  fe.bind(*element);
  return fe;
}

/**
 * \brief A test suite to test a set of bases.
 * \tparam numberOfBases Number of bases to be tested.
 * \tparam FEBases Type of different FEBases.
 * \param expectedNumberOfChildren Expected number of children for each basis.
 * \param expectedSize Expected size for each basis.
 * \param feBases A tuple of FEBases.
 * \param basisName Name of the types of bases.
 * \return Test suite object.
 */
template <size_t numberOfBases, typename... FEBases>
auto FEBaseAndIndicesTest(const std::array<int, numberOfBases>& expectedNumberOfChildren,
                          const std::array<int, numberOfBases>& expectedSize, const std::tuple<FEBases...>& feBases,
                          const std::string& basisName) {
  static_assert(std::tuple_size_v<std::tuple<FEBases...>> == numberOfBases,
                "Input size mismatch in FEBaseAndIndicesTest.");
  TestSuite t("FEBaseAndIndicesTest");
  Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<numberOfBases>()), [&](const auto i) {
    const auto febase         = std::get<i>(feBases);
    const std::string message = " for i = " + std::to_string(i) + " -> " + basisName;
    using MultiIndex          = typename std::remove_cvref_t<decltype(febase.localView())>::MultiIndex;
    static_assert(Dune::IsIndexable<MultiIndex>(), "MultiIndex must support operator[]");
    checkScalars(t, static_cast<int>(febase.localView().tree().degree()), expectedNumberOfChildren[i],
                 " Number of children is incorrect" + message);
    checkScalars(t, static_cast<int>(febase.size()), expectedSize[i], " Size is incorrect" + message);
    std::vector<MultiIndex> dofs;
    Ikarus::FEHelper::globalIndicesFromLocalView(febase.localView(), dofs);
    checkScalars(t, dofs.size(), febase.size(), " Size of dofs and basis is not equal" + message);
    t.check(!dofs.empty(), "dofs is empty" + message);
    std::ranges::sort(dofs);
    const bool hasDuplicates = std::adjacent_find(dofs.begin(), dofs.end()) == dofs.end();
    t.check(hasDuplicates) << "The sorted dofs vector has duplicates" + message;
  });
  return t;
}

template <bool useFlat>
auto FEBaseTest() {
  TestSuite t("FEBaseTest");
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox       = {4, 2};
  std::array<int, 2> elementsPerDirection = {1, 1};
  const auto grid                         = std::make_shared<Grid>(bbox, elementsPerDirection);

  const auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;

  const auto firstOrderLagrangePreBasis  = lagrange<1>();
  const auto secondOrderLagrangePreBasis = lagrange<2>();

  const auto firstOrderLagrangePowerPreBasis  = power<2>(firstOrderLagrangePreBasis);
  const auto secondOrderLagrangePowerPreBasis = power<5>(secondOrderLagrangePreBasis);

  const auto scalarScalarCompositePreBasis = composite(firstOrderLagrangePreBasis, secondOrderLagrangePreBasis);
  const auto scalarPowerCompositePreBasis  = composite(firstOrderLagrangePreBasis, secondOrderLagrangePowerPreBasis);

  const auto powerScalarCompositePreBasis = composite(firstOrderLagrangePowerPreBasis, secondOrderLagrangePreBasis);
  const auto powerPowerCompositePreBasis = composite(firstOrderLagrangePowerPreBasis, secondOrderLagrangePowerPreBasis);

  const auto scalarCompositePreBasis = composite(firstOrderLagrangePreBasis);
  const auto powerCompositePreBasis  = composite(firstOrderLagrangePowerPreBasis);

  const auto combinedPreBasis =
      composite(composite(scalarCompositePreBasis, powerPowerCompositePreBasis, firstOrderLagrangePreBasis),
                secondOrderLagrangePowerPreBasis);

  const auto compositePowerCombinedBasis = composite(power<10>(combinedPreBasis), secondOrderLagrangePreBasis);

  const auto lagrangeDGPowerPreBasis     = power<6>(lagrangeDG<3>());
  const auto nedelecScalarPreBasis       = nedelec<1, 1>();
  const auto raviartThomasScalarPreBasis = raviartThomas<2>();
  const auto specialCompositePreBasis =
      composite(lagrangeDGPowerPreBasis, nedelecScalarPreBasis, raviartThomasScalarPreBasis);

  /// Types of scalar bases
  const auto scalar1 = getFEBase<useFlat>(gridView, firstOrderLagrangePreBasis);
  const auto scalar2 = getFEBase<useFlat>(gridView, secondOrderLagrangePreBasis);

  /// Types of power bases
  const auto power1 = getFEBase<useFlat>(gridView, firstOrderLagrangePowerPreBasis);
  const auto power2 = getFEBase<useFlat>(gridView, secondOrderLagrangePowerPreBasis);

  /// Types of composite bases
  const auto composite1 = getFEBase<useFlat>(gridView, scalarScalarCompositePreBasis);
  const auto composite2 = getFEBase<useFlat>(gridView, scalarPowerCompositePreBasis);
  const auto composite3 = getFEBase<useFlat>(gridView, powerScalarCompositePreBasis);
  const auto composite4 = getFEBase<useFlat>(gridView, powerPowerCompositePreBasis);
  const auto composite5 = getFEBase<useFlat>(gridView, scalarCompositePreBasis);
  const auto composite6 = getFEBase<useFlat>(gridView, powerCompositePreBasis);
  const auto composite7 = getFEBase<useFlat>(gridView, combinedPreBasis);
  const auto composite8 = getFEBase<useFlat>(gridView, compositePowerCombinedBasis);

  /// Types of special bases
  const auto special1 = getFEBase<useFlat>(gridView, lagrangeDGPowerPreBasis);
  const auto special2 = getFEBase<useFlat>(gridView, nedelecScalarPreBasis);
  const auto special3 = getFEBase<useFlat>(gridView, raviartThomasScalarPreBasis);
  const auto special4 = getFEBase<useFlat>(gridView, specialCompositePreBasis);

  t.subTest(FEBaseAndIndicesTest<2>({0, 0}, {4, 9}, std::make_tuple(scalar1, scalar2), "Scalar basis"));
  t.subTest(FEBaseAndIndicesTest<2>({2, 5}, {8, 45}, std::make_tuple(power1, power2), "Power Basis"));
  t.subTest(FEBaseAndIndicesTest<6>(
      {2, 2, 2, 2, 1, 1}, {13, 49, 17, 53, 4, 8},
      std::make_tuple(composite1, composite2, composite3, composite4, composite5, composite6), "Composite Basis"));
  t.subTest(FEBaseAndIndicesTest<2>({2, 2}, {106, 1069}, std::make_tuple(composite7, composite8),
                                    "Combined Composite Basis"));
  t.subTest(FEBaseAndIndicesTest<1>({6}, {96}, std::make_tuple(special1), "Lagrange DG Basis"));
  t.subTest(FEBaseAndIndicesTest<1>({0}, {4}, std::make_tuple(special2), "Nedelec Basis"));
  t.subTest(FEBaseAndIndicesTest<1>({0}, {24}, std::make_tuple(special3), "Raviart Thomas Basis"));
  t.subTest(FEBaseAndIndicesTest<1>({3}, {124}, std::make_tuple(special4), "Special Composite Basis"));

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(FEBaseTest<true>());
  t.subTest(FEBaseTest<false>());
  return t.exit();
}
