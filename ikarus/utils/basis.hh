// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file basis.hh
 * \brief Wrapper around Dune-functions global basis

 */

#pragma once

#include <utility>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

#include <ikarus/utils/flatprebasis.hh>

namespace Ikarus {

/**
 * \brief Wrapper class for a hierarchical basis constructed from a pre-basis.
 *
 * This class provides a convenient wrapper for a hierarchical basis constructed from a pre-basis.
 * It contains both the original untouched basis and a flattened version of the basis.
 * \ingroup  utils
 * \tparam PB The type of the pre-basis used to construct the hierarchical basis.
 */
template <typename PB>
class BasisHandler
{
public:
  using PreBasis = PB;                          ///< The type of the untouched pre basis
  using GridView = typename PreBasis::GridView; ///< The type of the grid view
  using UntouchedBasis =
      decltype(Dune::Functions::DefaultGlobalBasis(std::declval<PreBasis>())); ///< The type of the untouched basis
  using FlatBasis = decltype(Dune::Functions::DefaultGlobalBasis(
      Ikarus::flatPreBasis(std::declval<PreBasis>()))); ///< The type of the flattened  basis

  /**
   * \brief Constructs a Basis object from a pre-basis.
   *
   * Constructs a Basis object from a given pre-basis, creating both the untouched and flat versions of the basis.
   *
   * \param pb The pre-basis used to construct the hierarchical basis.
   */
  explicit BasisHandler(const PreBasis& pb)
      : bb_{Dune::Functions::DefaultGlobalBasis(pb)},
        fb_{Dune::Functions::DefaultGlobalBasis(Ikarus::flatPreBasis(pb))} {}

  /**
   * \brief Returns a reference to the flat version of the basis.
   *
   * \return A reference to the flat version of the basis.
   */
  auto& flat() { return fb_; }

  /**
   * \brief Returns a reference to the untouched version of the basis.
   *
   * \return A reference to the untouched version of the basis.
   */
  auto& untouched() { return bb_; }

  /**
   * \brief Returns a const reference to the flat version of the basis.
   *
   * \return A const reference to the flat version of the basis.
   */
  const auto& flat() const { return fb_; }

  /**
   * \brief Returns a const reference to the untouched version of the basis.
   *
   * \return A const reference to the untouched version of the basis.
   */
  const auto& untouched() const { return bb_; }

  /**
   * \brief Returns a const reference to the grid view associated with the untouched basis.
   *
   * \return A const reference to the grid view associated with the untouched basis.
   */
  const auto& gridView() const { return bb_.gridView(); }

  /**
   * \brief Returns a reference to the grid view associated with the untouched basis.
   *
   * \return A reference to the grid view associated with the untouched basis.
   */
  auto& gridView() { return bb_.gridView(); }

private:
  UntouchedBasis bb_;
  FlatBasis fb_;
};

/**
 * \brief Factory function to create a Basis object.
 *
 * Factory function that creates a Basis object from a given grid view and pre-basis factory.
 *
 * \tparam GV The type of the grid view.
 * \tparam PBF The type of the pre-basis factory function.
 * \param gv The grid view.
 * \param pb The pre-basis factory function.
 * \return A Basis object.
 */
template <typename GV, typename PBF>
auto makeBasis(const GV& gv, const PBF& pb) {
  auto preBasis = pb(gv);
  return BasisHandler(preBasis);
}

/**
 * \brief Factory function to create a Basis object from a DefaultGlobalBasis.
 *
 *
 * \tparam PB The type of the pre-basis associated with the DefaultGlobalBasis.
 * \param gb The DefaultGlobalBasis.
 * \return A Basis object.
 */
template <typename PB>
auto makeBasis(const Dune::Functions::DefaultGlobalBasis<PB>& gb) {
  return BasisHandler(gb.preBasis());
}
} // namespace Ikarus
