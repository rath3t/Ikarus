// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file fetraits.hh
 * \brief FETraits template structure for finite element traits.
 */

#pragma once

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/concepts.hh>

namespace Ikarus {

/**
 * \brief Traits for handling finite elements.
 *
 * \tparam BH The basis handler type for the finite element.
 * \tparam useFlat A boolean indicating if the type of the underlying basis is of the flat or the untouched version.
 * \tparam useRef Boolean indicating whether to use Eigen::Ref for VectorType and MatrixType.
 */
template <typename BH, bool useRef = false, bool useFlat = true>
struct FETraits
{
  /** \brief Type of the basis of the finite element */
  using BasisHandler = BH;

  /** \brief A bool to indicate if the provided basishandler should hand out the flat basis */
  static constexpr bool useFlatBasis = useFlat;

  /** \brief Type of the flat basis */
  using FlatBasis = typename BasisHandler::FlatBasis;

  /** \brief Type of the untouched basis */
  using UntouchedBasis = typename BasisHandler::UntouchedBasis;

  /** \brief Type of the basis version*/
  using Basis = std::conditional_t<useFlat, FlatBasis, UntouchedBasis>;

  /** \brief Type of the local view */
  using LocalView = typename Basis::LocalView;

  /** \brief Type of the grid view */
  using GridView = typename Basis::GridView;

  /** \brief Type of the grid element */
  using Element = typename LocalView::Element;

  /** \brief Type of the element geometry */
  using Geometry = typename Element::Geometry;

  /** \brief Type of the global index */
  using GlobalIndex = typename LocalView::MultiIndex;

  /** \brief Type used for coordinates */
  using ctype = double;

  /** \brief Dimension of the world space */
  static constexpr int worlddim = Geometry::coorddimension;

  /** \brief Dimension of the geometry */
  static constexpr int mydim = Element::mydimension;

  /** \brief Dimension of the grid */
  static constexpr int dimension = Element::dimension;

  /** \brief Type of the  coordinate */
  using GlobalCoordinates = Eigen::Matrix<ctype, worlddim, 1>;

  /** \brief Type of the ParameterSpace coordinate */
  using ParameterSpaceType = Eigen::Matrix<ctype, mydim, 1>;

  /** \brief Bool indicating whether the raw Eigen types should be used or wrapped with Eigen::Ref<..>. (Needed for
   * Python bindings) */
  static constexpr bool useEigenRef = useRef;

  /** \brief Type of the vector passed to calculateVector */
  template <typename ScalarType = ctype>
  using VectorType =
      std::conditional_t<useEigenRef, Eigen::Ref<Eigen::VectorX<ScalarType>>, Eigen::VectorX<ScalarType>&>;

  /** \brief Type of the vector passed to updateState */
  template <typename ScalarType = ctype>
  using VectorTypeConst =
      std::conditional_t<useEigenRef, const Eigen::Ref<Eigen::VectorX<ScalarType>>, const Eigen::VectorX<ScalarType>&>;

  /** \brief Type of the matrix  passed to calculateMatrix */
  template <typename ScalarType = ctype>
  using MatrixType =
      std::conditional_t<useEigenRef, Eigen::Ref<Eigen::MatrixX<ScalarType>>, Eigen::MatrixX<ScalarType>&>;
};

} // namespace Ikarus
