// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file easvariants.hh
 * \brief Definition of the EAS variants.
 * \ingroup  mechanics
 */

#pragma once

#include <ikarus/finiteelements/mechanics/strainenhancements/linearandglstrains.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::EAS {
namespace Impl {
  template <typename GEO, int myDim, StrainTags ESS>
  struct EASType;
  template <typename GEO>
  struct EASType<GEO, 3, StrainTags::linear>
  {
    using type = std::variant<E0<GEO>, E9<GEO>, E21<GEO>>;
  };
  template <typename GEO>
  struct EASType<GEO, 2, StrainTags::linear>
  {
    using type = std::variant<E0<GEO>, E4<GEO>, E5<GEO>, E7<GEO>>;
  };
  template <typename GEO, int myDim>
  struct EASType<GEO, myDim, StrainTags::greenLagrangian>
  {
    using type = EASType<GEO, myDim, StrainTags::linear>::type;
  };
} // namespace Impl

/**
 * \brief Wrapper around the EAS variant, contains helper functions
 * \tparam GEO the Geometry type of the udderlying grid element
 */
template <StrainTags ES, typename GEO>
struct EASVariant
{
  static constexpr int myDim = GEO::mydimension;
  using Variant              = Impl::EASType<GEO, myDim, ES>::type;
  static_assert((myDim == 2) or (myDim == 3), "EAS variants are only available 2D and 3D elements.");

  /**
   * \brief Executes a function F on the variant, if `enhancedStrainSize` is not zero
   * \tparam F the Type of the function F
   * \param f the function, which takes the eas element as an argument
   */
  template <typename F>
  void operator()(F&& f) const {
    std::visit([&]<typename EAST>(const EAST& easFunction) { f(easFunction); }, var_);
  }

  /**
   * \brief A helper function to get the number of EAS parameters.
   * \return Number of EAS parameters.
   */
  auto numberOfEASParameters() const {
    return std::visit([]<typename EAST>(const EAST&) -> int { return EAST::enhancedStrainSize; }, var_);
  }

  /**
   * \brief A helper function to identify if the formulation is purely displacement-based, i.e., number of EAS
   * parameters is zero.
   * \return A bool indicating if the formulation is displacement-based or not.
   */
  bool isDisplacmentBased() const { return numberOfEASParameters() == 0; }

  void setEASType(int numberOfEASParameters) {
    numberOfEASParameters_ = numberOfEASParameters;
    if (geometry_)
      createEASType();
  }
  void bind(const GEO& geometry) {
    geometry_ = std::make_optional<GEO>(geometry);
    createEASType();
  }

private:
  void createEASType() {
    const std::string& errorMessage = "The given EAS parameters are not available for enhancing " + toString(ES) +
                                      " strain measure for the " + std::to_string(myDim) + "D case.";

    if constexpr (ES == StrainTags::linear or ES == StrainTags::greenLagrangian) {
      if (numberOfEASParameters_ == 0) {
        var_ = E0(geometry_.value());
        return;
      }
      if constexpr (myDim == 2) {
        switch (numberOfEASParameters_) {
          case 4:
            var_ = E4(geometry_.value());
            break;
          case 5:
            var_ = E5(geometry_.value());
            break;
          case 7:
            var_ = E7(geometry_.value());
            break;
          default:
            DUNE_THROW(Dune::NotImplemented, errorMessage);
        }
      } else if constexpr (myDim == 3) {
        switch (numberOfEASParameters_) {
          case 9:
            var_ = E9(geometry_.value());
            break;
          case 21:
            var_ = E21(geometry_.value());
            break;
          default:
            DUNE_THROW(Dune::NotImplemented, errorMessage);
        }
      }
    }
  }
  std::optional<GEO> geometry_;
  Variant var_;
  int numberOfEASParameters_;
};
} // namespace Ikarus::EAS
