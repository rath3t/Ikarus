// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file linearsolver.hh
 * \brief Type-erased linear solver with templated scalar type
 */

#pragma once
#include <memory>
#include <type_traits>
#include <variant>

#include <dune/common/exceptions.hh>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <ikarus/utils/makeenum.hh>
namespace Ikarus {

/**
 * \enum SolverTypeTag
 * \brief Enumeration representing different solver types.
 * \details The prefix s and d stand for sparse and dense solvers and the second prefix i and d stand for iterative or
 * direct solvers for the sparse case
 */
MAKE_ENUM(SolverTypeTag, none, si_ConjugateGradient, si_LeastSquaresConjugateGradient, si_BiCGSTAB, sd_SimplicialLLT,
          sd_SimplicialLDLT, sd_SparseLU, sd_SparseQR, sd_CholmodSupernodalLLT, sd_UmfPackLU, sd_SuperLU,
          d_PartialPivLU, d_FullPivLU, d_HouseholderQR, d_ColPivHouseholderQR, d_FullPivHouseholderQR,
          d_CompleteOrthogonalDecomposition, d_LLT, d_LDLT);

/**
 * \enum MatrixTypeTag
 * \brief Enumeration representing different matrix types (Dense or Sparse).
 */
enum class MatrixTypeTag
{
  Dense,
  Sparse
};

/**
 * \class LinearSolverTemplate
 * \brief A type-erased class which wraps most of the linear solvers available in Eigen.
 * \tparam ST The scalar type of the linear system (default: double).
 * \ingroup solvers
 */
template <typename ST = double>
class LinearSolverTemplate
{
public:
  using ScalarType       = ST;
  using SparseMatrixType = Eigen::SparseMatrix<ScalarType>;
  using DenseMatrixType  = Eigen::MatrixX<ScalarType>;

  LinearSolverTemplate() = default;

  /**
   * \brief Constructor for LinearSolverTemplate.
   * \param solverTypeTag The solver type tag representing the type of the linear solver.
   */
  explicit LinearSolverTemplate(const SolverTypeTag& solverTypeTag);

  /**
   * \brief Destructor for LinearSolverTemplate.
   */
  ~LinearSolverTemplate() = default;

  friend void swap(LinearSolverTemplate& lhs, LinearSolverTemplate& rhs) noexcept {
    std::swap(lhs.solverimpl_, rhs.solverimpl_);
    std::swap(lhs.solverTypeTag_, rhs.solverTypeTag_);
  }
  /**
   * \brief Copy assignment operator.
   * \param other The LinearSolverTemplate to copy.
   * \return A reference to the assigned LinearSolverTemplate.
   */
  LinearSolverTemplate& operator=(const LinearSolverTemplate& other) {
    LinearSolverTemplate tmp(other);
    swap(*this, tmp);
    return *this;
  }

  /**
   * \brief Copy constructor.
   * \param other The LinearSolverTemplate to copy.
   */
  LinearSolverTemplate(const LinearSolverTemplate& other) { *this = LinearSolverTemplate(other.solverTypeTag_); }
  /**
   * \brief Move constructor.
   * \param other The LinearSolverTemplate to move.
   */
  LinearSolverTemplate(LinearSolverTemplate&& other) noexcept = default;
  /**
   * \brief Move assignment operator.
   * \param other The LinearSolverTemplate to move.
   * \return A reference to the assigned LinearSolverTemplate.
   */
  LinearSolverTemplate& operator=(LinearSolverTemplate&& other) noexcept = default;

private:
  struct SolverBase
  {
    virtual ~SolverBase() = default;
    virtual void analyzePattern(const DenseMatrixType&) const {}
    virtual void analyzePattern(const SparseMatrixType&)                                         = 0;
    virtual void factorize(const DenseMatrixType&)                                               = 0;
    virtual void factorize(const SparseMatrixType&)                                              = 0;
    virtual void compute(const SparseMatrixType&)                                                = 0;
    virtual void compute(const DenseMatrixType&)                                                 = 0;
    virtual void solve(Eigen::VectorX<ScalarType>& x, const Eigen::VectorX<ScalarType>&) const   = 0;
    virtual void solve(Eigen::MatrixX2<ScalarType>& x, const Eigen::MatrixX2<ScalarType>&) const = 0;
    virtual void solve(Eigen::MatrixX3<ScalarType>& x, const Eigen::MatrixX3<ScalarType>&) const = 0;
    virtual void solve(Eigen::MatrixX<ScalarType>& x, const Eigen::MatrixX<ScalarType>&) const   = 0;
  };

  template <typename Solver>
  struct SolverImpl : public SolverBase
  {
    using SolverBase::analyzePattern; // forward use of analyzePattern(DenseMatrixType)

    void analyzePattern(const SparseMatrixType& A) final {
      if constexpr (requires(Solver sol) { sol.analyzePattern(A); })
        solver.analyzePattern(A);
    }

    void factorize(const SparseMatrixType& A) final {
      if constexpr (requires(Solver sol) { sol.factorize(A); })
        solver.factorize(A);
    }

    // Dense Solvers do not have a factorize method therefore
    // our interface we just call compute for dense matrices
    void factorize(const DenseMatrixType& A) final {
      if constexpr (requires(Solver sol) { sol.compute(A); } && std::is_base_of_v<Eigen::SolverBase<Solver>, Solver>)
        solver.compute(A);
    }
    void compute(const SparseMatrixType& A) final {
      if constexpr (std::is_base_of_v<Eigen::SparseSolverBase<Solver>, Solver>)
        solver.compute(A);
      else
        DUNE_THROW(Dune::NotImplemented, "This solver does not support solving with sparse matrices.");
    }

    void compute(const DenseMatrixType& A) final {
      if constexpr (std::is_base_of_v<Eigen::SolverBase<Solver>, Solver>)
        solver.compute(A);
      else
        DUNE_THROW(Dune::NotImplemented, "This solver does not support solving with dense matrices.");
    }

    void solve(Eigen::VectorX<ScalarType>& x, const Eigen::VectorX<ScalarType>& b) const final { x = solver.solve(b); }

    void solve(Eigen::MatrixX2<ScalarType>& x, const Eigen::MatrixX2<ScalarType>& b) const final {
      x = solver.solve(b);
    }

    void solve(Eigen::MatrixX3<ScalarType>& x, const Eigen::MatrixX3<ScalarType>& b) const final {
      x = solver.solve(b);
    }

    void solve(Eigen::MatrixX<ScalarType>& x, const Eigen::MatrixX<ScalarType>& b) const final { x = solver.solve(b); }

    Solver solver;
  };

  std::unique_ptr<SolverBase> solverimpl_;
  SolverTypeTag solverTypeTag_{SolverTypeTag::none};

public:
  /**
   * \brief Compute the factorization of the matrix.
   * \tparam MatrixType The type of the matrix (DenseMatrixType or SparseMatrixType).
   * \param A The matrix for factorization.
   * \return A reference to the LinearSolverTemplate.
   */
  template <typename MatrixType>
  requires std::is_same_v<MatrixType, DenseMatrixType> || std::is_same_v<MatrixType, SparseMatrixType>
  inline LinearSolverTemplate& compute(const MatrixType& A) {
    solverimpl_->compute(A);
    return *this;
  }

  /**
   * \brief Analyze the pattern of the matrix.
   * \tparam MatrixType The type of the matrix (DenseMatrixType or SparseMatrixType).
   * \param A The matrix for pattern analysis.
   */
  template <typename MatrixType>
  requires std::is_same_v<MatrixType, DenseMatrixType> || std::is_same_v<MatrixType, SparseMatrixType>
  inline void analyzePattern(const MatrixType& A) {
    solverimpl_->analyzePattern(A);
  }

  /**
   * \brief Factorize the matrix.
   * \tparam MatrixType The type of the matrix (DenseMatrixType or SparseMatrixType).
   * \param A The matrix for factorization.
   */
  template <typename MatrixType>
  requires std::is_same_v<MatrixType, DenseMatrixType> || std::is_same_v<MatrixType, SparseMatrixType>
  inline void factorize(const MatrixType& A) {
    solverimpl_->factorize(A);
  }

  /**
   * \brief Solve the linear system for a vector.
   * \param x The solution vector.
   * \param b The right-hand side vector.
   */
  void solve(Eigen::VectorX<ScalarType>& x, const Eigen::VectorX<ScalarType>& b) { solverimpl_->solve(x, b); }

  /**
   * \brief Solve the linear system for a `n` times `3` matrix.
   * \param x The solution matrix.
   * \param b The right-hand side matrix.
   */
  void solve(Eigen::MatrixX3<ScalarType>& x, const Eigen::MatrixX3<ScalarType>& b) { solverimpl_->solve(x, b); }
  /**
   * \brief Solve the linear system for a `n` times `2` matrix.
   * \param x The solution matrix.
   * \param b The right-hand side matrix.
   */
  void solve(Eigen::MatrixX2<ScalarType>& x, const Eigen::MatrixX2<ScalarType>& b) { solverimpl_->solve(x, b); }

  /**
   * \brief Solve the linear system for a `n` times `n` matrix.
   * \param x The solution matrix.
   * \param b The right-hand side matrix.
   */
  void solve(Eigen::MatrixX<ScalarType>& x, const Eigen::MatrixX<ScalarType>& b) { solverimpl_->solve(x, b); }
};

using LinearSolver = LinearSolverTemplate<double>;
extern template class LinearSolverTemplate<double>;
} // namespace Ikarus
