// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file febase.hh
 * \brief Contains the FE class, which is used as a base class for all finite elements.
 * It provides information about the local view of a finite element.
 * It also contains the PreFE class which acts as a convenient wrapper
 * for different type traits needed by the FE class.
 * \ingroup finiteelements
 */

#pragma once

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/fetraits.hh>
#include <ikarus/finiteelements/mixin.hh>

namespace Ikarus {

template <typename PreFE, template <typename, typename> class... Skills>
class FE;

/**
 * \brief PreFE struct acts as a convenient wrapper for the FE class to access different type traits.
 *
 * \tparam BH The type of the basis handler.
 * \tparam useFlat A boolean indicating if the type of the underlying basis is of the flat or the untouched version.
 * \tparam useEigenRef A boolean flag indicating whether to use Eigen references.
 */
template <typename BH, bool useFlat = true, bool useEigenRef = false>
struct PreFE
{
  using BasisHandler                      = BH;
  static constexpr int worldDim           = BasisHandler::Basis::worlddim;
  static constexpr bool useEigenReference = useEigenRef;
  static constexpr bool useFlatBasis      = useFlat;
  using Traits                            = FETraits<BH, useEigenRef, useFlat>;

  template <template <typename, typename> class... Skills>
  using FE = FE<PreFE, Skills...>;
};

namespace Impl {
  // Since base classes are initialized, in declaration order before member variables, we have to make sure the
  // localview_ object of FE is defined. To do this we add this artificial inheritance by inheriting below first from
  // FEInit, which initializes the localview object first, siche we the other FEMixin base class constructors are called
  template <typename PreFE, typename FE>
  struct FEInit
  {
    using Traits       = PreFE::Traits;              ///< Type traits
    using BasisHandler = Traits::BasisHandler;       ///< Type of the basisHandler.
    using LocalView    = typename Traits::LocalView; ///< Type of the local view.
    friend FE;
    FEInit(const BasisHandler& basisHandler)
        : localView_{[&]() -> LocalView {
            if constexpr (Traits::useFlatBasis)
              return basisHandler.flat().localView();
            else
              return basisHandler.untouched().localView();
          }()} {}

  private:
    LocalView localView_;
  };
} // namespace Impl

/**
 * \brief FE class is a base class for all finite elements.
 *
 * \details A single class which can be used by the finite elements to get different information
 * from the local view of the element. It works for both the flat and untouched version of the basis.
 *
 * \tparam PreFE The type of the  pre finite element.
 * \tparam Skills Variadic template for skill arguments possessed by the FE.
 */
template <typename PreFE, template <typename, typename> class... Skills>
class FE : private Impl::FEInit<PreFE, FE<PreFE, Skills...>>, public FEMixin<PreFE, Skills...>
{
  friend Impl::FEInit<PreFE, FE>;

protected:
  using Mixin = FEMixin<PreFE, Skills...>; ///< Type of the FE mixin.
  friend Mixin;

public:
  using Traits       = PreFE::Traits;                ///< Type traits
  using BasisHandler = Traits::BasisHandler;         ///< Type of the basisHandler.
  using LocalView    = typename Traits::LocalView;   ///< Type of the local view.
  using GridView     = typename Traits::GridView;    ///< Type of the global view.
  using GlobalIndex  = typename Traits::GlobalIndex; ///< Type of the global index.
  using GridElement  = typename Traits::Element;     ///< Type of the grid element.
  using typename Mixin::Requirement;

  static constexpr int myDim    = Traits::mydim;
  static constexpr int worldDim = Traits::worlddim;
  using PreTuple                = std::tuple<typename Skills<PreFE, FE>::Pre...>;

  /**
   * \brief Constructor for the FE class.
   *
   * \param basisHandler The basis handler.
   * \param skillsArgs Skill arguments.
   */
  explicit FE(const BasisHandler& basisHandler, typename Skills<PreFE, FE>::Pre... skillsArgs)
      : Impl::FEInit<PreFE, FE>(basisHandler),
        FEMixin<PreFE, Skills...>(std::forward<typename Skills<PreFE, FE>::Pre>(skillsArgs)...) {}

  /**
   * \brief Convenient function to bind the local view to the element.
   * @param  element The element to be bounded
   */
  void bind(const GridElement& element) {
    this->localView_.bind(element);
    Mixin::bind();
  }

  /**
   * \brief Get the size of the local view.
   * \return The size of the local view.
   */
  [[nodiscard]] constexpr size_t size() const { return this->localView_.size(); }

  /**
   * \brief Get the grid element associated with the local view.
   * \return The grid element.
   */
  const GridElement& gridElement() const { return this->localView_.element(); }

  /**
   * \brief Get the const reference to the local view.
   * \return The const reference to the local view.
   */
  const LocalView& localView() const { return this->localView_; }

  /**
   * \brief Get the reference to the local view.
   * \return The reference to the local view.
   */
  LocalView& localView() { return this->localView_; }
};
} // namespace Ikarus
