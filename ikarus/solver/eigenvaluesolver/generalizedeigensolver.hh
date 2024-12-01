// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file generalizedeigensolver.hh
 * \brief Implementation of wrapper classes for solving a generalized eigenvalue problem
 */

#pragma once

#include <optional>
#include <utility>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymGEigsSolver.h>

#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/makeenum.hh>

namespace Ikarus {

/**
 * \brief Construct a new make enum object
 */
MAKE_ENUM(EigenValueSolverType, Spectra, Eigen);

template <EigenValueSolverType SolverType, Concepts::DenseOrSparseEigenMatrix MT>
struct GeneralizedSymEigenSolver
{
};

/**
 * \brief This class implements a wrapper to the Spectra generalized eigen solver for real symmetric matrices, i.e. to
 * solve \f$ Ax = \lambda Bx\f$, where A is symmetric and B is positive definite. It calculates the full spectrum of
 * eigenvalues.
 * \details Under the hood it uses the Spectra::SymGEigsSolver with Cholesky decomposition for B. As this class can only
 * compute up to \f$ n-1 \f$ smallest or greatest eigenvalues, we use two different solvers for the n/2 smallest and n/2
 * greatest eigenvalues. The matrices are shared throughout the solvers so no extra copy is being made.
 *
 * \tparam MT the used Matrix Type, can be a sparse or dense Eigen::Matrix
 */
template <Concepts::DenseOrSparseEigenMatrix MT>
struct GeneralizedSymEigenSolver<EigenValueSolverType::Spectra, MT>
{
  using ScalarType              = typename MT::Scalar;
  static constexpr bool isDense = Concepts::EigenMatrix<MT>;

  using MatrixType = MT;
  using ProductType =
      std::conditional_t<isDense, Spectra::DenseSymMatProd<ScalarType>, Spectra::SparseSymMatProd<ScalarType>>;
  using CholeskyType =
      std::conditional_t<isDense, Spectra::DenseCholesky<ScalarType>, Spectra::SparseCholesky<ScalarType>>;

  using SolverType = Spectra::SymGEigsSolver<ProductType, CholeskyType, Spectra::GEigsMode::Cholesky>;

  /**
   * \brief Construct a new General Sym Eigen Solver object.
   *
   * \tparam MATA deduced type of passed matrix A.
   * \tparam MATB deduced type of passed matrix B.
   * \param A matrix A.
   * \param B matrix B.
   */
  template <typename MATA, typename MATB>
  requires(Concepts::DenseOrSparseEigenMatrix<std::remove_cvref_t<MATA>>)
  GeneralizedSymEigenSolver(MATA&& A, MATB&& B)
      : nev_(A.rows()),
        nevsPartition_(static_cast<Eigen::Index>(std::ceil(static_cast<double>(nev_) / 2)),
                       static_cast<Eigen::Index>(nev_ - std::ceil(static_cast<double>(nev_) / 2))),
        aOP_(std::forward<MATA>(A)),
        bOP_(std::forward<MATB>(B)),
        solverSmallest_(aOP_, bOP_, nevsPartition_.first, std::min(nevsPartition_.second * 2, nev_)),
        solverGreatest_(aOP_, bOP_, nevsPartition_.second, nevsPartition_.second * 2) {
    if ((A.cols() != B.cols()) or (A.rows() != B.rows()) or (A.cols() != B.cols()))
      DUNE_THROW(Dune::IOError, "GeneralizedSymEigenSolver: The passed matrices should have the same size");
    eigenvalues_.resize(nev_);
    eigenvectors_.resize(A.rows(), nev_);
    assert(nevsPartition_.first + nevsPartition_.second == nev_);
  }

  /**
   * \brief Construct a new General Sym Eigen Solver object.
   *
   * \tparam AssemblerA the type of the assembler for matrix A.
   * \tparam AssemblerB the type of the assembler for matrix B.
   * \param assemblerA assembler for matrix A.
   * \param assemblerB assembler for matrix B.
   */
  template <Concepts::FlatAssembler AssemblerA, Concepts::FlatAssembler AssemblerB>
  GeneralizedSymEigenSolver(const std::shared_ptr<AssemblerA>& assemblerA,
                            const std::shared_ptr<AssemblerB>& assemblerB)
      : GeneralizedSymEigenSolver(assemblerA->matrix(), assemblerB->matrix()) {}

  /**
   * \brief Starts the computation of the eigenvalue solver
   *
   * \param tolerance given tolerance for iterative eigenvalue solving (default: 1e-10)
   * \param maxit givenn maximum iterations for eigenvalue solving (default: 1000)
   * \return true solving was successful
   * \return false solving was not successful
   */
  bool compute(ScalarType tolerance = 1e-10, Eigen::Index maxit = 1000) {
    solverSmallest_.init();
    solverSmallest_.compute(Spectra::SortRule::SmallestAlge, 1000, 1e-10, Spectra::SortRule::SmallestAlge);

    eigenvalues_.head(nevsPartition_.first)      = solverSmallest_.eigenvalues();
    eigenvectors_.leftCols(nevsPartition_.first) = solverSmallest_.eigenvectors();

    solverGreatest_.init();
    solverGreatest_.compute(Spectra::SortRule::LargestAlge, 1000, 1e-10, Spectra::SortRule::SmallestAlge);

    eigenvalues_.tail(nevsPartition_.second)       = solverGreatest_.eigenvalues();
    eigenvectors_.rightCols(nevsPartition_.second) = solverGreatest_.eigenvectors();

    computed_ = solverSmallest_.info() == Spectra::CompInfo::Successful and
                solverGreatest_.info() == Spectra::CompInfo::Successful;

    return computed_;
  }

  /**
   * \brief Returns the eigenvalues of the generalized eigenvalue problem
   *
   * \return Eigen::VectorXd vector of eigenvalues
   */
  auto& eigenvalues() const {
    assertCompute();
    return eigenvalues_;
  }

  /**
   * \brief Returns the eigenvectors of the generalized eigenvalue problem
   *
   * \return auto matrix with the eigevectors as columns
   */
  auto& eigenvectors() const {
    assertCompute();
    return eigenvectors_;
  }

  /** \brief Returns the number of eigenvalues of the problem */
  Eigen::Index nev() const { return nev_; }

private:
  Eigen::Index nev_;
  std::pair<Eigen::Index, Eigen::Index> nevsPartition_;
  ProductType aOP_;
  CholeskyType bOP_;
  SolverType solverSmallest_;
  SolverType solverGreatest_;

  Eigen::VectorXd eigenvalues_;
  Eigen::MatrixXd eigenvectors_;

  bool computed_{};

  void assertCompute() const {
    if (not computed_)
      DUNE_THROW(Dune::IOError, "Eigenvalues and -vectors not yet computed, please call compute() first");
  }
};

/**
 * \brief This class implements a wrapper to the Eigen generalized eigen solver for real symmetric matrices, i.e. to
 * solve \f$ Ax = \lambda Bx\f$, where A is symmetric and B is positive definite. A and B have to be dense matrices.
 * \details Under the hood it uses the Eigen::GeneralizedSelfAdjointEigenSolver
 *
 * \tparam MT the used Matrix Type, can be a dense Eigen::Matrix
 */
template <Concepts::EigenMatrix MT>
struct GeneralizedSymEigenSolver<EigenValueSolverType::Eigen, MT>
{
  using ScalarType = typename MT::Scalar;
  using MatrixType = MT;
  using SolverType = Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType>;

  template <typename MATA, typename MATB>
  requires(Concepts::EigenMatrix<std::remove_cvref_t<MATA>>)
  GeneralizedSymEigenSolver(MATA&& A, MATB&& B)
      : matA_(std::forward<MATA>(A)),
        matB_(std::forward<MATB>(B)),
        solver_(A.size()) {
    if ((A.cols() != B.cols()) or (A.rows() != B.rows()) or (A.cols() != B.cols()))
      DUNE_THROW(Dune::IOError, "GeneralizedSymEigenSolver: The passed matrices should have the same size");
  }

  template <Concepts::FlatAssembler AssemblerA, Concepts::FlatAssembler AssemblerB>
  GeneralizedSymEigenSolver(const std::shared_ptr<AssemblerA>& assemblerA, const std::shared_ptr<AssemblerB>& assemblerB)
      : GeneralizedSymEigenSolver(assemblerA->matrix(), assemblerB->matrix()) {}

  /**
   * \brief Starts the computation of the eigenvalue solver
   *
   * \param options defaults to Eigen::ComputeEigenvectors, can be set to Eigen::EigenvaluesOnly, accessing eigenvectors
   * in that case will results in an error
   * \return true solving was successful
   * \return false solving was not successful
   */
  bool compute(int options = Eigen::ComputeEigenvectors) {
    solver_.compute(matA_, matB_, options);
    computed_ = true;
    return true; // SelfAdjointEigenSolver will always be successful if prerequisites are met
  }

  /**
   * \brief Returns the eigenvalues of the generalized eigenvalue problem
   *
   * \return Reference to the vector of eigenvalues
   */
  auto& eigenvalues() const {
    assertCompute();
    return solver_.eigenvalues();
  }

  /**
   * \brief Returns the eigenvectors of the generalized eigenvalue problem
   *
   * \param _nev optionally specify how many eigenvectors are requested
   * \return Reference to the matrix with the eigevectors as columns
   */
  auto& eigenvectors() const {
    assertCompute();
    return solver_.eigenvectors();
  }

private:
  MatrixType matA_;
  MatrixType matB_;
  SolverType solver_;
  bool computed_{};

  void assertCompute() const {
    if (not computed_)
      DUNE_THROW(Dune::IOError, "Eigenvalues and -vectors not yet computed, please call compute() first");
  }
};

template <EigenValueSolverType tag, Concepts::FlatAssembler AS1, Concepts::FlatAssembler AS2>
requires(std::same_as<typename AS1::MatrixType, typename AS2::MatrixType> &&
         not(tag == EigenValueSolverType::Eigen && Concepts::SparseEigenMatrix<typename AS1::MatrixType>))
auto makeGeneralizedSymEigenSolver(const std::shared_ptr<AS1>& as1, const std::shared_ptr<AS2> as2) {
  using MatrixType        = typename AS1::MatrixType;
  constexpr auto isSparse = Concepts::SparseEigenMatrix<MatrixType>;
  using SolverType        = std::conditional_t<isSparse, GeneralizedSymEigenSolver<tag, MatrixType>,
                                               GeneralizedSymEigenSolver<tag, MatrixType>>;

  return SolverType{as1, as2};
}

/**
 * \brief This class implements a wrapper to the Spectra generalized eigen solver for real symmetric matrices, i.e. to
 * solve \f$ Ax = \lambda Bx\f$, where A is symmetric and B is positive definite. It calculates a selection of
 * eigenvalues. At most \f$ n - 1 \f$, where \f$ n \f$ is the rows/cols of the matrices.
 * \details Under the hood it uses the Spectra::SymGEigsSolver with Cholesky decomposition for B
 *
 * \tparam MT the used Matrix Type, can be a sparse or dense Eigen::Matrix
 */
template <Concepts::DenseOrSparseEigenMatrix MT>
struct PartialGeneralizedSymEigenSolver
{
  using ScalarType              = typename MT::Scalar;
  static constexpr bool isDense = Concepts::EigenMatrix<MT>;

  using MatrixType = MT;
  using ProductType =
      std::conditional_t<isDense, Spectra::DenseSymMatProd<ScalarType>, Spectra::SparseSymMatProd<ScalarType>>;
  using CholeskyType =
      std::conditional_t<isDense, Spectra::DenseCholesky<ScalarType>, Spectra::SparseCholesky<ScalarType>>;

  using SolverType = Spectra::SymGEigsSolver<ProductType, CholeskyType, Spectra::GEigsMode::Cholesky>;

  template <typename MATA, typename MATB>
  requires(Concepts::DenseOrSparseEigenMatrix<std::remove_cvref_t<MATA>>)
  PartialGeneralizedSymEigenSolver(MATA&& A, MATB&& B, Eigen::Index nev)
      : nev_(nev),
        aOP_(std::forward<MATA>(A)),
        bOP_(std::forward<MATB>(B)),
        solver_(aOP_, bOP_, nev, 2 * nev <= A.cols() ? 2 * nev : A.cols()) {
    if ((A.cols() != B.cols()) or (A.rows() != B.rows()) or (A.cols() != B.cols()))
      DUNE_THROW(Dune::IOError, "PartialGeneralizedSymEigenSolver: The passed matrices should have the same size");
  }

  template <Concepts::FlatAssembler AssemblerA, Concepts::FlatAssembler AssemblerB>
  PartialGeneralizedSymEigenSolver(const std::shared_ptr<AssemblerA>& assemblerA,
                                   const std::shared_ptr<AssemblerB>& assemblerB, Eigen::Index nev)
      : PartialGeneralizedSymEigenSolver(assemblerA->matrix(), assemblerB->matrix(), nev) {}

  /**
   * \brief Starts the computation of the eigenvalue solver
   *
   * \param selection defines which n eigenvalues should be returned (i.e. smallest or greatesst).
   * \param sortRule defines how these eigenvalues should be sorted.
   * \param tolerance given tolerance for iterative eigenvalue solving (default: 1e-10).
   * \param maxit givenn maximum iterations for eigenvalue solving (default: 1000).
   * \return true solving was successful.
   * \return false solving was not successful.
   */
  bool compute(Spectra::SortRule selection = Spectra::SortRule::SmallestAlge,
               Spectra::SortRule sortRule = Spectra::SortRule::SmallestAlge, ScalarType tolerance = 1e-10,
               Eigen::Index maxit = 1000) {
    solver_.init();
    solver_.compute(selection, tolerance, maxit, sortRule);

    computed_ = solver_.info() == Spectra::CompInfo::Successful;
    return computed_;
  }

  /**
   * \brief Returns the eigenvalues of the generalized eigenvalue problem
   *
   * \return Eigen::VectorXd vector of eigenvalues
   */
  Eigen::VectorXd eigenvalues() const {
    assertCompute();
    return solver_.eigenvalues();
  }

  /**
   * \brief Returns the eigenvectors of the generalized eigenvalue problem
   *
   * \param _nev optionally specify how many eigenvectors are requested
   * \return auto matrix with the eigevectors as columns
   */
  auto eigenvectors(std::optional<Eigen::Index> _nev = std::nullopt) const {
    assertCompute();
    return solver_.eigenvectors(_nev.value_or(nev_));
  }

  Eigen::Index nev() const { return nev_; }

private:
  Eigen::Index nev_;
  ProductType aOP_;
  CholeskyType bOP_;
  SolverType solver_;

  bool computed_{};

  void assertCompute() const {
    if (not computed_)
      DUNE_THROW(Dune::IOError, "Eigenvalues and -vectors not yet computed, please call compute() first");
  }
};

template <Concepts::FlatAssembler AS1, Concepts::FlatAssembler AS2>
auto makePartialGeneralizedSymEigenSolver(const std::shared_ptr<AS1>& as1, const std::shared_ptr<AS2> as2) {
  using MatrixType = typename AS1::MatrixType;
  using ScalarType = typename AS1::MatrixType::Scalar;
  using SolverType = PartialGeneralizedSymEigenSolver<MatrixType>;

  return SolverType{as1, as2};
}

#ifndef DOXYGEN
template <Concepts::FlatAssembler AS1, Concepts::FlatAssembler AS2>
PartialGeneralizedSymEigenSolver(std::shared_ptr<AS1> as1, std::shared_ptr<AS2> as2,
                                 int nev) -> PartialGeneralizedSymEigenSolver<typename AS1::MatrixType>;
#endif
} // namespace Ikarus
