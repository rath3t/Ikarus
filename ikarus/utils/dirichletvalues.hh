// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file dirichletvalues.hh
 * \brief Definition of DirichletValues class for handling Dirichlet boundary conditions.
 *
 */

#pragma once
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <functional>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

#include <Eigen/Core>

#include <autodiff/forward/real/real.hpp>

#include <ikarus/utils/concepts.hh>

namespace Dune {
template <class K, int N>
class FieldVector;
}

namespace Ikarus {

/**
 * \brief A helper struct to derive the SizeType of the underlying container.
 */
template <class>
struct DeriveSizeType;

template <Concepts::EigenVector T>
struct DeriveSizeType<T>
{
  using SizeType = Eigen::Index;
};

template <>
struct DeriveSizeType<std::vector<bool>>
{
  using SizeType = std::vector<bool>::size_type;
};

/**
 * \brief Class for handling Dirichlet boundary conditions in Ikarus.
 * \ingroup  utils
 * \details The DirichletValues class provides functionalities for fixing degrees of freedom and storing inhomogeneous
 * Dirichlet boundary conditions. It supports fixing degrees of freedom using various callback functions and
 * stores functions for inhomogeneous Dirichlet boundary conditions.
 *
 * \tparam B Type of the finite element basis
 * \tparam FC Type for storing Dirichlet flags (default is std::vector<bool>)
 */
template <typename B, typename FC = std::vector<bool>>
class DirichletValues
{
public:
  using Basis                         = std::remove_cvref_t<B>;
  using FlagsType                     = FC;
  static constexpr int worldDimension = Basis::GridView::dimensionworld;
  using BackendType                   = decltype(Dune::Functions::istlVectorBackend(std::declval<FlagsType&>()));
  using SizeType                      = typename DeriveSizeType<FlagsType>::SizeType;
  explicit DirichletValues(const B& basis)
      : basis_{basis},
        dirichletFlagsBackend_{dirichletFlags_} {
    dirichletFlagsBackend_.resize(basis_);
    std::fill(dirichletFlags_.begin(), dirichletFlags_.end(), false);
    inhomogeneousBoundaryVectorDummy_.setZero(static_cast<Eigen::Index>(basis_.size()));
  }

  /**
   * \brief Function to fix (set boolean values to true or false) degrees of freedom on the boundary.
   *
   * This function takes a callback function, which will be called with the boolean vector of fixed boundary
   * degrees of freedom and the usual arguments of `Dune::Functions::forEachBoundaryDOF`.
   *
   * \param f A callback function
   * \param treePath An optional argument specifying a tree path to a subspacebasis, e.g. Dune::Indices::_0
   */
  template <typename F, typename TreePath = Dune::TypeTree::HybridTreePath<>>
  void fixBoundaryDOFs(F&& f, TreePath&& treePath = {}) {
    using namespace Dune::Functions;
    using SubSpaceLocalView =
        typename std::remove_cvref_t<decltype(subspaceBasis(basis_, std::forward<TreePath>(treePath)))>::LocalView;

    if constexpr (Concepts::IsFunctorWithArgs<F, BackendType, typename Basis::MultiIndex>) {
      auto lambda = [&](auto&& indexGlobal) { f(dirichletFlagsBackend_, indexGlobal); };
      Dune::Functions::forEachBoundaryDOF(subspaceBasis(basis_, std::forward<TreePath>(treePath)), lambda);
    } else if constexpr (Concepts::IsFunctorWithArgs<F, BackendType, int, SubSpaceLocalView>) {
      auto lambda = [&](auto&& localIndex, auto&& localView) { f(dirichletFlagsBackend_, localIndex, localView); };
      Dune::Functions::forEachBoundaryDOF(subspaceBasis(basis_, std::forward<TreePath>(treePath)), lambda);
    } else if constexpr (Concepts::IsFunctorWithArgs<F, BackendType, int, SubSpaceLocalView,
                                                     typename Basis::GridView::Intersection>) {
      auto lambda = [&](auto&& localIndex, auto&& localView, auto&& intersection) {
        f(dirichletFlagsBackend_, localIndex, localView, intersection);
      };
      Dune::Functions::forEachBoundaryDOF(subspaceBasis(basis_, std::forward<TreePath>(treePath)), lambda);
    } else {
      static_assert(Dune::AlwaysFalse<F>(), "fixBoundaryDOFs: A function with this signature is not supported");
    }
  }

  /**
   * \brief Function to fix (set boolean values to true or false) degrees of freedom.
   *
   * This function takes a callback function, which will be called with the stored function space basis and
   * the boolean vector of fixed boundary degrees of freedom.
   *
   * \param f A callback function
   */
  template <typename F>
  void fixDOFs(F&& f) {
    f(basis_, dirichletFlagsBackend_);
  }

  /**
   * \brief Fixes and unfixes (set boolean value to true or false) a specific degree of freedom
   *
   * \param i An index indicating the DOF number to be fixed
   * \param flag Boolean indicating whether the DOF should fixed or not
   */
  template <typename MultiIndex>
  requires(not std::integral<MultiIndex>)
  void setSingleDOF(const MultiIndex i, bool flag) {
    dirichletFlagsBackend_[i] = flag;
  }

  /**
   * \brief Fixes or unfixes (set boolean value to true or false) a specific degree of freedom
   *
   * \param i An index indicating the DOF number to be fixed
   * \param flag Boolean indicating whether the DOF should fixed or not
   */

  void setSingleDOF(std::size_t i, bool flag)
  requires(std::same_as<typename Basis::MultiIndex, Dune::Functions::FlatMultiIndex<size_t>>)
  {
    dirichletFlags_[i] = flag;
  }

  /**
   * \brief Resets all degrees of freedom
   */
  void reset() {
    std::fill(dirichletFlags_.begin(), dirichletFlags_.end(), false);
    inhomogeneousBoundaryVectorDummy_.setZero(static_cast<Eigen::Index>(basis_.size()));
  }

  /* \brief Returns the local basis object */
  const auto& basis() const { return basis_; }

  /* \brief Returns a boolean values, if the give multiIndex is constrained */
  template <typename MultiIndex>
  requires(not std::integral<MultiIndex>)
  [[nodiscard]] bool isConstrained(const MultiIndex& multiIndex) const {
    return dirichletFlagsBackend_[multiIndex];
  }

  /* \brief Returns a boolean values, if the i-th degree of freedom is constrained */
  [[nodiscard]] bool isConstrained(std::size_t i) const
  requires(std::same_as<typename Basis::MultiIndex, Dune::Functions::FlatMultiIndex<size_t>>)
  {
    return dirichletFlags_[i];
  }

  /* \brief Returns how many degrees of freedoms are fixed */
  auto fixedDOFsize() const { return std::ranges::count(dirichletFlags_, true); }

  /* \brief Returns how many degrees of freedom there are */
  auto size() const { return dirichletFlags_.size(); }

  /* \brief Returns the underlying dof flag container */
  auto& container() const { return dirichletFlags_; }

  /**
   * \brief Function to insert a function for inhomogeneous Dirichlet boundary conditions.
   *
   * \details This function takes a callback function, which will be called with the current coordinate vector
   * and the scalar load factor. It creates internally the first derivative of the passed function and
   * stores them simultaneously.
   *
   * \param f A callback function
   */
  template <typename F>
  void storeInhomogeneousBoundaryCondition(F&& f) {
    auto derivativeLambda = [&](const auto& globalCoord, const double& lambda) {
      autodiff::real lambdaDual = lambda;
      lambdaDual[1]             = 1; // Setting the derivative in lambda direction to 1
      return derivative(f(globalCoord, lambdaDual));
    };
    dirichletFunctions_.push_back({f, derivativeLambda});
  }

  /**
   * \brief Function to evaluate all stored inhomogeneous Dirichlet boundary functions.
   *
   * This function evaluates all stored inhomogeneous Dirichlet boundary functions at all positions where the
   * corresponding degrees of freedom are true.
   *
   * \param xIh The vector where the interpolated result should be stored
   * \param lambda The load factor
   */
  void evaluateInhomogeneousBoundaryCondition(Eigen::VectorXd& xIh, const double& lambda) {
    inhomogeneousBoundaryVectorDummy_.setZero();
    xIh.resizeLike(inhomogeneousBoundaryVectorDummy_);
    xIh.setZero();
    for (auto& f : dirichletFunctions_) {
      interpolate(
          basis_, inhomogeneousBoundaryVectorDummy_,
          [&](const auto& globalCoord) { return f.value(globalCoord, lambda); }, dirichletFlagsBackend_);
      xIh += inhomogeneousBoundaryVectorDummy_;
    }
  }

  /**
   * \brief Function to evaluate all stored inhomogeneous Dirichlet boundary derivative functions.
   *
   * This function evaluates all stored inhomogeneous Dirichlet boundary derivative functions at all positions
   * where the corresponding degrees of freedom are true.
   *
   * \param xIh The vector where the interpolated result should be stored
   * \param lambda The load factor
   */
  void evaluateInhomogeneousBoundaryConditionDerivative(Eigen::VectorXd& xIh, const double& lambda) {
    inhomogeneousBoundaryVectorDummy_.setZero();
    xIh.resizeLike(inhomogeneousBoundaryVectorDummy_);
    xIh.setZero();
    for (auto& f : dirichletFunctions_) {
      interpolate(
          basis_, inhomogeneousBoundaryVectorDummy_,
          [&](const auto& globalCoord) { return f.derivative(globalCoord, lambda); }, dirichletFlagsBackend_);
      xIh += inhomogeneousBoundaryVectorDummy_;
    }
  }

private:
  Eigen::VectorXd inhomogeneousBoundaryVectorDummy_;
  Basis basis_;
  FlagsType dirichletFlags_;
  BackendType dirichletFlagsBackend_;
  struct DirichletFunctions
  {
    using Signature = std::function<Eigen::Vector<double, worldDimension>(
        const Dune::FieldVector<double, worldDimension>&, const double&)>;
    Signature value;
    Signature derivative;
  };
  std::vector<DirichletFunctions> dirichletFunctions_;
};

} // namespace Ikarus
