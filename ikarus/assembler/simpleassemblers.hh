// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file simpleassemblers.hh
 * \brief Defines several assemblers for finite element assembly
 */

#pragma once

#include <utility>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <ikarus/assembler/interface.hh>
#include <ikarus/finiteelements/ferequirements.hh>

namespace Ikarus {
/**
 * \class ScalarFlatAssembler
 * \brief ScalarFlatAssembler assembles scalar quantities.
 * \ingroup assembler
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 */
template <typename FEC, typename DV>
class ScalarFlatAssembler : public ScalarAssembler<ScalarFlatAssembler<FEC, DV>, FEC, DV, double>,
                            public FlatAssemblerBase<FEC, DV>
{
protected:
  using Base = FlatAssemblerBase<FEC, DV>; ///< Type alias for the base class.
  friend ScalarAssembler<ScalarFlatAssembler, FEC, DV, double>;

  using Base::Base;

public:
  using typename Base::Basis;
  using typename Base::DirichletValuesType;
  using typename Base::FEContainer;
  using typename Base::FERequirement;
  using typename Base::GlobalIndex;
  using typename ScalarAssembler<ScalarFlatAssembler, FEC, DV, double>::ScalarType;

protected:
  ScalarType& getScalarImpl(const FERequirement& feRequirements, ScalarAffordance affordance);

  ScalarType scal_{0.0};
};

#ifndef DOXYGEN
template <class T, class DirichletValuesType>
ScalarFlatAssembler(T&& fes, const DirichletValuesType& dirichletValues) -> ScalarFlatAssembler<T, DirichletValuesType>;
#endif

/**
 * \class VectorFlatAssembler
 * \brief VectorFlatAssembler assembles vector quantities using a flat basis Indexing strategy.
 *\ingroup assembler
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 */
template <typename FEC, typename DV>
class VectorFlatAssembler : public ScalarFlatAssembler<FEC, DV>,
                            public VectorAssembler<VectorFlatAssembler<FEC, DV>, FEC, DV, Eigen::VectorXd>
{
protected:
  using Base = ScalarFlatAssembler<FEC, DV>; ///< Type alias for the base class.
  friend VectorAssembler<VectorFlatAssembler, FEC, DV, Eigen::VectorXd>;

  using Base::Base;

public:
  using typename Base::Basis;
  using typename Base::DirichletValuesType;
  using typename Base::FEContainer;
  using typename Base::FERequirement;
  using typename Base::GlobalIndex;

  using typename Base::ScalarType;
  using typename VectorAssembler<VectorFlatAssembler, FEC, DV, Eigen::VectorXd>::VectorType;

protected:
  void assembleRawVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance, VectorType& assemblyVec);
  VectorType& getRawVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance);
  VectorType& getVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance);

  VectorType& getReducedVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance);

  VectorType vecRaw_{}; ///< Raw vector without changes for dirichlet degrees of freedom
  VectorType vec_{};    ///< Vector quantity.
  VectorType vecRed_{}; ///< Reduced vector quantity.
};

#ifndef DOXYGEN
template <class T, class DirichletValuesType>
VectorFlatAssembler(T&& fes, const DirichletValuesType& dirichletValues) -> VectorFlatAssembler<T, DirichletValuesType>;
#endif

/**
 * \class SparseFlatAssembler
 * \brief SparseFlatAssembler assembles matrix quantities using a flat basis Indexing strategy.
 * The matrix is stored in a sparse matrix format. This format is exploited during the assembly process.
 * \ingroup assembler
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 */
template <typename FEC, typename DV>
class SparseFlatAssembler : public MatrixAssembler<SparseFlatAssembler<FEC, DV>, FEC, DV, Eigen::SparseMatrix<double>>,
                            public VectorFlatAssembler<FEC, DV>
{
protected:
  using Base = VectorFlatAssembler<FEC, DV>; ///< Type alias for the base class.
  friend MatrixAssembler<SparseFlatAssembler, FEC, DV, Eigen::SparseMatrix<double>>;

  using Base::Base;

public:
  using typename Base::Basis;
  using typename Base::DirichletValuesType;
  using typename Base::FEContainer;
  using typename Base::FERequirement;
  using typename Base::GlobalIndex;
  using typename Base::ScalarType;
  using typename Base::VectorType;
  using typename MatrixAssembler<SparseFlatAssembler, FEC, DV, Eigen::SparseMatrix<double>>::MatrixType;

private:
  void assembleRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance, MatrixType& assemblyMat);

protected:
  MatrixType& getRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);
  MatrixType& getMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);
  MatrixType& getReducedMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);

private:
  /** Calculates the non-zero entries in the full sparse matrix and passes them to the underlying Eigen sparse matrix.
   */
  void createOccupationPattern(MatrixType& assemblyMat);

  /** Calculates the non-zero entries in the sparse matrix and passes them to the underlying Eigen sparse matrix.
   * The size of the matrix has the size of the free degrees of freedom. */
  void createReducedOccupationPattern(MatrixType& assemblyMat);

  /** Saves the degree of freedom indices of each element in the vector elementLinearIndices. */
  void createLinearDOFsPerElement(MatrixType& assemblyMat);

  /** Saves the degree of freedom indices of each element in the vector elementLinearIndices but excludes fixed
   * degrees of freedom. */
  void createLinearDOFsPerElementReduced(MatrixType& assemblyMat);

  /** Pre-processes the raw sparse matrix before assembly. */
  void preProcessSparseMatrix(MatrixType& assemblyMat);

  /** Pre-processes the reduced sparse matrix before assembly. */
  void preProcessSparseMatrixReduced(MatrixType& assemblyMat);

  MatrixType spMatRaw_;     ///< Raw sparse matrix without changes for dirichlet degrees of freedom.
  MatrixType spMat_;        ///< Sparse matrix.
  MatrixType spMatReduced_; ///< Reduced sparse matrix.
  std::vector<std::vector<Eigen::Index>>
      elementLinearIndices_; ///< Vector storing indices of matrix entries in linear storage
  std::vector<std::vector<Eigen::Index>>
      elementLinearReducedIndices_; ///< Vector storing indices of matrix entries in linear storage
  bool sparsePreProcessorRaw_{}, sparsePreProcessor_{},
      sparsePreProcessorReduced_{}; ///< flags that store if the sparsity pattern construction happened
};

#ifndef DOXYGEN
template <class FEC, class DV>
SparseFlatAssembler(FEC&& fes, const DV& dirichletValues) -> SparseFlatAssembler<FEC, DV>;
#endif

template <typename FEC, typename DV>
auto makeSparseFlatAssembler(FEC&& fes, const DV& dirichletValues) {
  return std::make_shared<SparseFlatAssembler<FEC, DV>>(std::forward<FEC>(fes), dirichletValues);
}

/**
 * \class DenseFlatAssembler
 * \brief DenseFlatAssembler assembles matrix quantities using a flat basis Indexing strategy.
 * The matrix is stored in a dense matrix format. This format is exploited during the assembly process.
 * \ingroup assembler
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 * \note Requires Ikarus::Concepts::FlatIndexBasis<BasisEmbedded>.
 */
template <typename FEC, typename DV>
class DenseFlatAssembler : public MatrixAssembler<DenseFlatAssembler<FEC, DV>, FEC, DV, Eigen::MatrixXd>,
                           public VectorFlatAssembler<FEC, DV>
{
protected:
  using Base = VectorFlatAssembler<FEC, DV>; ///< Type alias for the base class.
  friend MatrixAssembler<DenseFlatAssembler, FEC, DV, Eigen::MatrixXd>;

  using Base::Base;

public:
  using typename Base::Basis;
  using typename Base::DirichletValuesType;
  using typename Base::FEContainer;
  using typename Base::FERequirement;
  using typename Base::GlobalIndex;
  using typename Base::ScalarType;
  using typename Base::VectorType;
  using typename MatrixAssembler<DenseFlatAssembler, FEC, DV, Eigen::MatrixXd>::MatrixType;

private:
  void assembleRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance, MatrixType& assemblyMat);

protected:
  MatrixType& getRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);
  MatrixType& getMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);
  MatrixType& getReducedMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);

private:
  MatrixType matRaw_{}; ///< Raw dense matrix for assembly.
  MatrixType mat_{};    ///< Dense matrix quantity.
  MatrixType matRed_{}; ///< Reduced dense matrix quantity.
};

#ifndef DOXYGEN
// https://en.cppreference.com/w/cpp/language/class_template_argument_deduction
template <class FEC, class DV>
DenseFlatAssembler(FEC&& fes, const DV& dirichletValues) -> DenseFlatAssembler<FEC, DV>;
#endif

template <typename FEC, typename DV>
auto makeDenseFlatAssembler(FEC&& fes, const DV& dirichletValues) {
  return std::make_shared<DenseFlatAssembler<FEC, DV>>(std::forward<FEC>(fes), dirichletValues);
}
} // namespace Ikarus

#include "simpleassemblers.inl"
