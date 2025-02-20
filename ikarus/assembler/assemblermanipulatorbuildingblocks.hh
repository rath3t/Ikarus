// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file assemblerwrappers.hh
 * \brief Defines the base classes for scalar, vector, and matrix wrappers to assemblers
 */

#pragma once

#include <utility>

#include <dune/functions/backends/istlvectorbackend.hh>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/assembler/interface.hh>
#include <ikarus/finiteelements/feconfiguration.hh>
#include <ikarus/finiteelements/fehelper.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/dirichletvalues.hh>

namespace Ikarus {

/**
 * \brief Base class for a wrapper to a scalar assembler.
 * \tparam Wrapper Type of the wrapper to an assembler
 * \tparam Assembler Type of the assembler
 */
template <typename Wrapper, typename Assembler>
struct ScalarManipulator
{
  using WrappedAssembler = Wrapper;
  using FEC              = typename Assembler::FEContainer;
  using DV               = typename Assembler::DirichletValuesType;
  using FEConfiguration  = typename Assembler::FEConfiguration;
  using ScalarType       = typename Assembler::ScalarType;
  using Interface        = ScalarAssembler<WrappedAssembler, FEC, DV, ScalarType>;
  friend Interface;

  using FunctionType = std::function<void(const Assembler&, const FEConfiguration&, ScalarAffordance, ScalarType&)>;

  /**
   * \brief A helper function to add functions that can be used to manipulate the assembled quantity.
   * \tparam F Type of the function
   * \param f A function that manipulates the assembled quantity.
   */
  template <typename F>
  requires Concepts::IsFunctorWithArgs<F, const Assembler&, const FEConfiguration&, ScalarAffordance, ScalarType&>
  void bind(F&& f) {
    sfs.emplace_back(std::forward<F>(f));
  }
  std::vector<FunctionType> sfs;

protected:
  ScalarType& getScalarImpl(const FEConfiguration& feconfiguration, ScalarAffordance affordance) {
    ScalarType& sca = underlying().base_getScalarImpl(feconfiguration, affordance);
    for (const auto& sf : sfs)
      sf(underlying().base(), feconfiguration, affordance, sca);
    return sca;
  }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const WrappedAssembler&>(*this); }
  auto& underlying() { return static_cast<WrappedAssembler&>(*this); }
};

/**
 * \brief Base class for a wrapper to a vector assembler.
 * \tparam Wrapper Type of the wrapper to an assembler
 * \tparam Assembler Type of the assembler
 */
template <typename Wrapper, typename Assembler>
struct VectorManipulator

{
  using WrappedAssembler = Wrapper;
  using FEC              = typename Assembler::FEContainer;
  using DV               = typename Assembler::DirichletValuesType;
  using FEConfiguration  = typename Assembler::FEConfiguration;

  using VectorType = typename Assembler::VectorType;
  using Interface  = VectorAssembler<WrappedAssembler, typename Assembler::FEContainer,
                                     typename Assembler::DirichletValuesType, typename Assembler::VectorType>;
  friend Interface;
  using FunctionType =
      std::function<void(const Assembler&, const FEConfiguration&, VectorAffordance, DBCOption, VectorType&)>;

  /**
   * \brief A helper function to add functions that can be used to manipulate the assembled quantity.
   * \tparam F Type of the function
   * \param f A function that manipulates the assembled quantity.
   */
  template <typename F>
  requires Concepts::IsFunctorWithArgs<F, const Assembler&, const FEConfiguration&, VectorAffordance, DBCOption,
                                       VectorType&>
  void bind(F&& f) {
    vfs.emplace_back(std::forward<F>(f));
  }
  std::vector<FunctionType> vfs;

protected:
  VectorType& getRawVectorImpl(const FEConfiguration& feconfiguration, VectorAffordance affordance) {
    VectorType& vec = underlying().base_getRawVectorImpl(feconfiguration, affordance);
    for (const auto& vf : vfs)
      vf(underlying().base(), feconfiguration, affordance, DBCOption::Raw, vec);
    return vec;
  }

  VectorType& getVectorImpl(const FEConfiguration& feconfiguration, VectorAffordance affordance) {
    VectorType& vec = underlying().base_getVectorImpl(feconfiguration, affordance);
    for (const auto& vf : vfs)
      vf(underlying().base(), feconfiguration, affordance, DBCOption::Full, vec);
    return vec;
  }

  VectorType& getReducedVectorImpl(const FEConfiguration& feconfiguration, VectorAffordance affordance) {
    VectorType& vec = underlying().base_getReducedVectorImpl(feconfiguration, affordance);
    for (const auto& vf : vfs)
      vf(underlying().base(), feconfiguration, affordance, DBCOption::Reduced, vec);
    return vec;
  }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const WrappedAssembler&>(*this); }
  auto& underlying() { return static_cast<WrappedAssembler&>(*this); }
};

/**
 * \brief Base class for a wrapper to a matrix assembler.
 * \tparam Wrapper Type of the wrapper to an assembler
 * \tparam Assembler Type of the assembler
 */
template <typename Wrapper, typename Assembler>
struct MatrixManipulator
{
  using WrappedAssembler = Wrapper;
  using FEC              = typename Assembler::FEContainer;
  using DV               = typename Assembler::DirichletValuesType;
  using FEConfiguration  = typename Assembler::FEConfiguration;

  using MatrixType = typename Assembler::MatrixType;
  using Interface  = MatrixAssembler<WrappedAssembler, typename Assembler::FEContainer,
                                     typename Assembler::DirichletValuesType, typename Assembler::MatrixType>;
  friend Interface;
  using FunctionType =
      std::function<void(const Assembler&, const FEConfiguration&, MatrixAffordance, DBCOption, MatrixType&)>;

  /**
   * \brief A helper function to add functions that can be used to manipulate the assembled quantity.
   * \tparam F Type of the function
   * \param f A function that manipulates the assembled quantity.
   */
  template <typename F>
  requires Concepts::IsFunctorWithArgs<F, const Assembler&, const FEConfiguration&, MatrixAffordance, DBCOption,
                                       MatrixType&>
  void bind(F&& f) {
    mfs.emplace_back(std::forward<F>(f));
  }
  std::vector<FunctionType> mfs;

protected:
  MatrixType& getRawMatrixImpl(const FEConfiguration& feconfiguration, MatrixAffordance affordance) {
    MatrixType& mat = underlying().base_getRawMatrixImpl(feconfiguration, affordance);
    for (const auto& mf : mfs)
      mf(underlying().base(), feconfiguration, affordance, DBCOption::Raw, mat);
    return mat;
  }

  MatrixType& getMatrixImpl(const FEConfiguration& feconfiguration, MatrixAffordance affordance) {
    MatrixType& mat = underlying().base_getMatrixImpl(feconfiguration, affordance);
    for (const auto& mf : mfs)
      mf(underlying().base(), feconfiguration, affordance, DBCOption::Full, mat);
    return mat;
  }

  MatrixType& getReducedMatrixImpl(const FEConfiguration& feconfiguration, MatrixAffordance affordance) {
    MatrixType& mat = underlying().base_getReducedMatrixImpl(feconfiguration, affordance);
    for (const auto& mf : mfs)
      mf(underlying().base(), feconfiguration, affordance, DBCOption::Reduced, mat);
    return mat;
  }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const WrappedAssembler&>(*this); }
  auto& underlying() { return static_cast<WrappedAssembler&>(*this); }
};
} // namespace Ikarus
