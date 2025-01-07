// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file easvariants.hh
 * \brief Definition of the EAS variants.
 * \ingroup  mechanics
 */

#pragma once

#include <ikarus/finiteelements/mechanics/eas/eas2d.hh>
#include <ikarus/finiteelements/mechanics/eas/eas3d.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::EAS {
namespace Impl {
  /**
   * \brief EASVariants structure holding variants of different Enhanced Assumed Strains (EAS).
   *
   * \details EASVariants holds the different variants of EAS formulations for
   * both 2D and 3D elements.
   *
   * \tparam GEO The geometry type.
   */
  template <typename GEO>
  struct Variants
  {
    static constexpr int worldDim = GEO::coorddimension;
    static_assert((worldDim == 2) or (worldDim == 3), "EAS variants are only available 2D and 3D elements.");

    /** \brief Variant type holding different EAS formulations for 2D elements. */
    using EAS2D = std::variant<E0<GEO>, Q1E4<GEO>, Q1E5<GEO>, Q1E7<GEO>>;

    /** \brief Variant type holding different EAS formulations for 3D elements. */
    using EAS3D = std::variant<E0<GEO>, H1E9<GEO>, H1E21<GEO>>;

    /** \brief Type of the EAS variant depending on the grid dimension. */
    using type = std::conditional_t<worldDim == 2, EAS2D, EAS3D>;
  };

  /**
   * \brief Wrapper around the EAS variant, contains helper functions
   * \tparam Geometry the Geometry type of the udderlying grid element
   */
  template <typename Geometry>
  struct EASVariant
  {
    using Variant = Impl::Variants<Geometry>::type;

    /**
     * \brief Executes a function F on the variant, if `enhancedStrainSize` is not zero
     * \tparam F the Type of the function F
     * \param f the function, which takes the eas element as an argument
     */
    template <typename F>
    void operator()(F&& f) const {
      std::visit(
          [&]<typename EAST>(const EAST& easFunction) {
            if constexpr (EAST::enhancedStrainSize != 0)
              f(easFunction);
          },
          var_);
    }

    auto numberOfEASParameters() const {
      return std::visit([]<typename EAST>(const EAST&) -> int { return EAST::enhancedStrainSize; }, var_);
    }
    bool isDisplacmentBased() const { return numberOfEASParameters() == 0; }
    void setEASType(int numberOfEASParameters) {
      numberOfEASParameters_ = numberOfEASParameters;
      if (geometry_)
        createEASType();
    }
    void bind(const Geometry& geometry) {
      geometry_ = std::make_optional<Geometry>(geometry);
      createEASType();
    }

  private:
    void createEASType() {
      if (numberOfEASParameters_ == 0) {
        var_ = E0(geometry_.value());
        return;
      }

      if constexpr (Geometry::mydimension == 2) {
        switch (numberOfEASParameters_) {
          case 4:
            var_ = Q1E4(geometry_.value());
            break;
          case 5:
            var_ = Q1E5(geometry_.value());
            break;
          case 7:
            var_ = Q1E7(geometry_.value());
            break;
          default:
            DUNE_THROW(Dune::NotImplemented, "The given EAS parameters are not available for the 2D case.");
        }
      } else if constexpr (Geometry::mydimension == 3) {
        switch (numberOfEASParameters_) {
          case 9:
            var_ = H1E9(geometry_.value());
            break;
          case 21:
            var_ = H1E21(geometry_.value());
            break;
          default:
            DUNE_THROW(Dune::NotImplemented, "The given EAS parameters are not available for the 3D case.");
        }
      }
    }
    std::optional<Geometry> geometry_;
    Variant var_;
    int numberOfEASParameters_;
  };
} // namespace Impl
} // namespace Ikarus::EAS
