//
// Created by Alex on 26.05.2021.
//

#pragma once

#include <ranges>
#include <unordered_set>

#include <dune/common/hybridutilities.hh>
#include <dune/geometry/dimension.hh>

#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Grids/GridData.h>
#include <ikarus/Variables/InterfaceVariable.h>
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/Variables/VariableVector.h>
#include <ikarus/utils/utils/algorithms.h>
#include <ikarus/utils/utils/hashs.h>

namespace Ikarus::FEManager {

  template <typename FEContainer, typename GridViewType>
  requires Concepts::HasgetEntityVariableTuple<typename std::decay_t<FEContainer>::value_type> || Concepts::
      HasFreegetEntityVariableTuple<typename std::decay_t<FEContainer>::value_type>
  class DefaultFEManager {
  public:
    using GridDataType = GridData<typename GridViewType::IndexSetType>;
    DefaultFEManager(FEContainer& feContainer, GridViewType& gv,
                     std::optional<std::reference_wrapper<GridDataType>> gridData = std::nullopt)
        : gridView_{&gv}, gridData_{gridData}, varVec{feContainer} {}

    size_t numberOfDegreesOfFreedom() { return varVec.correctionSize(); }

    auto elementDofs() { return varVec.elementDofs(); };
    void addData(GridDataType& gridData) { gridData_ = gridData; }

    auto elementVariables() { return varVec.elementVariables(); };

    auto elementDofVectorSize() { return varVec.elementDofVectorSize(); };

    auto elementIndicesVariableTuple() {
      return varVec.transform_viewOverElements([&](auto& fe) {
        auto elementsVariables = varVec.variablesOfSingleElement(fe);
        auto elementsIndices   = varVec.dofIndicesOfSingleElement(fe);
        return std::make_tuple(std::ref(fe), elementsIndices, elementsVariables);
      });
    }

    auto elementIndicesVariableDataTuple() {
      if (gridData_.has_value())
        return varVec.transform_viewOverElements([&](auto& fe) {
          auto elementsVariables = varVec.variablesOfSingleElement(fe);
          auto elementsIndices   = varVec.dofIndicesOfSingleElement(fe);
          auto elementData       = gridData_.value().get().getAllSubEntityDataOfFE(fe);

          return std::make_tuple(std::ref(fe), elementsIndices, elementsVariables, elementData);
        });
      else
        throw std::logic_error(
            "GridData is Empty this method should not be called. Call elementDofsVariableTuple() instead");
    }

    auto& getVariables() { return varVec; };
    auto getGridView() { return gridView_; };

  private:
    GridViewType* gridView_;
    std::optional<std::reference_wrapper<GridDataType>> gridData_;
    Ikarus::Variable::VariableVector<FEContainer> varVec;
  };
}  // namespace Ikarus::FEManager