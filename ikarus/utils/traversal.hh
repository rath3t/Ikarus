// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file traversal.hh
 * \brief Contains functions to traverse through a tree to its different nodes
 */

#pragma once
#include <ikarus/utils/concepts.hh>

namespace Ikarus::utils {

/**
 * \brief A function which loops over all the nodes of a tree and performs different actions
 * for a power node (with leaf node as child) and a leaf node depending on the corresponding functor passed.
 *
 * \details This function is inspired from the function Dune::TypeTree::Detail::forEachNode
 * available in dune/typetree/traversal.hh
 *
 * \tparam T Type of tree.
 * \tparam TreePath Type of the HybridTreePath.
 * \tparam PowerFunc Type of the functor called for a power node.
 * \tparam LeafFunc Type of the functor called for a leaf node.
 * \param tree Tree whose nodes are to be looped over.
 * \param treePath The tree path which handles index values.
 * \param powerFunc A functor to be called for a power node.
 * \param leafFunc A functor to be called for a leaf node.
 */
template <class T, class TreePath, class PowerFunc, class LeafFunc>
void forEachLeafOrPowerLeafNode(T&& tree, TreePath&& treePath, PowerFunc&& powerFunc, LeafFunc&& leafFunc) {
  using Tree = std::decay_t<T>;
  if constexpr (Tree::isLeaf) {
    leafFunc(tree, treePath);
  } else if constexpr (Tree::isPower) {
    if constexpr (Tree::template Child<Dune::Indices::_0>::Type::isLeaf) {
      powerFunc(tree, treePath);
    } else {
      for (std::size_t i = 0; i < tree.degree(); ++i) {
        auto childTreePath = Dune::TypeTree::push_back(treePath, i);
        forEachLeafOrPowerLeafNode(tree.child(i), childTreePath, powerFunc, leafFunc);
      }
    }
  } else {
    auto indices = std::make_index_sequence<Tree::degree()>{};
    Dune::Hybrid::forEach(indices, [&](auto i) {
      auto childTreePath = Dune::TypeTree::push_back(treePath, i);
      forEachLeafOrPowerLeafNode(tree.child(i), childTreePath, powerFunc, leafFunc);
    });
  }
}

/**
 * \brief A helper function that helps in traversing over the local coordinates of an element and
 * call a user-desired function
 * \see Dune book page 314 for details
 * \tparam size Size of the global nodal coordinate vector
 * \tparam LV Type of the local view
 * \tparam F Type of the functor that traverses over the local coordinate of an element
 * \param localView Local view bounded to an element
 * \param f A function that traverses over the local coordinate of an element
 */
template <typename LV, typename F, int size = LV::Element::Geometry::coorddimension>
requires(std::convertible_to<F, std::function<bool(int, Dune::FieldVector<double, size> &&)>>)
void forEachLagrangeNodePosition(const LV& localView, F&& f) {
  static_assert(Concepts::LagrangeNode<std::remove_cvref_t<decltype(localView.tree().child(0))>>,
                "forEachLagrangeNodePosition is only supported for Lagrange power basis");
  assert(localView.bound() && "The local view must be bound to an element");
  const auto& localFE = localView.tree().child(0).finiteElement();
  std::vector<Dune::FieldVector<double, size>> lagrangeNodeCoords;
  lagrangeNodeCoords.resize(localFE.size());
  std::vector<double> out;
  for (int i = 0; i < size; i++) {
    auto ithCoord = [&i](const Dune::FieldVector<double, size>& x) { return x[i]; };
    localFE.localInterpolation().interpolate(ithCoord, out);
    for (std::size_t j = 0; j < out.size(); j++)
      lagrangeNodeCoords[j][i] = out[j];
  }
  for (int nodeNumber = 0; auto& nCoord : lagrangeNodeCoords)
    if (f(nodeNumber++, std::move(nCoord)))
      break;
}
} // namespace Ikarus::utils
