//
// Created by Alex on 19.05.2021.
//

#pragma once

#include <concepts>
#include <numeric>

#include "dune/common/classname.hh"

namespace Ikarus::Concepts {

  /**
   * \brief Concept of a variable
   *
   *
   */
  template <typename VariableType>
  concept Variable = requires(VariableType var, VariableType var2, typename VariableType::CorrectionType correction,
                              typename VariableType::CoordinateType value, int i) {
    typename VariableType::ctype;
    VariableType::valueSize;
    VariableType::correctionSize;
    VariableType::tagvalue;
    typename VariableType::CoordinateType;
    typename VariableType::CorrectionType;
    { var.getValue() } -> std::same_as<typename VariableType::CoordinateType>;
    { var.setValue(value) } -> std::same_as<void>;
    { var[i] } -> std::same_as<typename VariableType::ctype&>;
    { var.update(correction) } -> std::same_as<void>;
    { var == var2 } -> std::same_as<bool>;
    var.getTag();
  };

}  // namespace Ikarus::Concepts
