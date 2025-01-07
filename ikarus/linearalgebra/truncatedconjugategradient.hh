// Original File: https://gitlab.com/libeigen/eigen/-/blob/master/Eigen/src/IterativeLinearSolvers/ConjugateGradient.h
// SPDX-FileCopyrightText: 2011-2014 Copyright (C) Gael Guennebaud <gael.guennebaud@inria.fr>
// SPDX-License-Identifier: MPL-2.0
// Modifications:
// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file truncatedconjugategradient.hh
 * \brief Definition of TruncatedConjugateGradient class for solving linear systems using truncated conjugate gradient
 * method.
 *
 * This file defines the TruncatedConjugateGradient class, which is an iterative solver for solving linear systems
 * using the truncated conjugate gradient method. It includes modifications to the original Eigen library's
 * ConjugateGradient.h file.
 *
 */

#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Eigen {
enum class TCGStopReason : int
{
  negativeCurvature,
  exceededTrustRegion,
  reachedTargetResidualKappaLinear,
  reachedTargetResidualThetaSuperLinear,
  maximumInnerIterations,
  modelIncreased
};
template <typename Scalar>
struct TCGInfo
{
  TCGStopReason stop_tCG    = TCGStopReason::maximumInnerIterations;
  Scalar Delta              = 100000;
  Scalar kappa              = 0.1;
  Scalar theta              = 1.0;
  Eigen::Index mininner     = 1;
  Eigen::Index maxinner     = 1000;
  Eigen::Index numInnerIter = 0;

  void initRuntimeOptions(int _num_dof_solve) {
    maxinner = _num_dof_solve * 2;
    Delta    = 100000; // typical distance of manifold
    // Delta0 = Delta_bar/8.0;
  }
};
namespace internal {

  /**
   * @internal
   * \brief Low-level truncated conjugate gradient algorithm.
   * \tparam MatrixType Type of the matrix A.
   * \tparam Rhs Type of the right-hand side vector b.
   * \tparam Dest Type of the solution vector x.
   * \tparam Preconditioner Type of the preconditioner.
   * \param mat The matrix A.
   * \param rhs The right-hand side vector b.
   * \param x On input and initial solution, on output the computed solution.
   * \param precond A preconditioner being able to efficiently solve for an approximation of Ax=b (regardless of b).
   * \param iters On input the max number of iterations, on output the number of performed iterations.
   * \param tol_error On input the tolerance error, on output an estimation of the relative error.
   * \param _info Information about the truncated conjugate gradient algorithm.
   */
  template <typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
  void truncated_conjugate_gradient(const MatrixType& mat, const Rhs& rhs, Dest& x, const Preconditioner& precond,
                                    Eigen::Index& iters, typename Dest::RealScalar& tol_error,
                                    TCGInfo<typename Dest::RealScalar>& _info) {
    using std::abs;
    using std::sqrt;
    typedef typename Dest::RealScalar RealScalar;
    typedef typename Dest::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, 1> VectorType;

    RealScalar tol = tol_error;
    Index maxIters = iters;

    Index n = mat.cols();

    VectorType residual = rhs - mat * x; // initial residual

    RealScalar rhsNorm2             = rhs.norm();
    const RealScalar considerAsZero = (std::numeric_limits<RealScalar>::min)();

    if (rhsNorm2 <= considerAsZero) {
      x.setZero();
      iters     = 0;
      tol_error = 0;
      return;
    }
    RealScalar threshold     = numext::maxi(tol * tol * rhsNorm2 * rhsNorm2, considerAsZero);
    RealScalar residualNorm2 = residual.norm();
    if (residualNorm2 * residualNorm2 < threshold) {
      iters     = 0;
      tol_error = (residualNorm2 / rhsNorm2);
      return;
    }

    double e_Pd     = 0.0;
    double e_Pe_new = 0.0;
    double e_Pe     = x.squaredNorm();
    double d_Pd;
    double d_Hd;
    VectorType p(n);
    p = precond.solve(residual); // initial search direction
                                 // bool coutflag=true;
    _info.stop_tCG = TCGStopReason::maximumInnerIterations;
    VectorType z(n), tmp(n);
    RealScalar absNew = numext::real(residual.dot(p)); // the square of the absolute value of r scaled by invM
    d_Pd              = absNew;
    Index i           = 1;
    while (i < maxIters) {
      tmp.noalias() = mat * p; // the bottleneck of the algorithm
      d_Hd          = p.dot(tmp);
      Scalar alpha  = absNew / d_Hd; // the amount we travel on dir

      e_Pe_new = e_Pe + 2.0 * alpha * e_Pd + alpha * alpha * d_Pd;

      if (d_Hd <= 0 || e_Pe_new >= _info.Delta * _info.Delta) // negative curvature or execdet trustregion
      {
        double tau = (-e_Pd + sqrt(e_Pd * e_Pd + d_Pd * (_info.Delta * _info.Delta - e_Pe))) / d_Pd;

        x += tau * p;
        if (d_Hd <= 0)
          _info.stop_tCG = TCGStopReason::negativeCurvature;
        else
          _info.stop_tCG = TCGStopReason::exceededTrustRegion;

        break;
      }
      e_Pe = e_Pe_new;
      x += alpha * p;          // update solution
      residual -= alpha * tmp; // update residual

      residualNorm2 = residual.norm();

      if (i >= _info.mininner &&
          residualNorm2 <= rhsNorm2 * std::min(rhsNorm2, _info.kappa)) // missing pow(rhsNorm2,_info.theta
      {
        // Residual is small enough to quit
        if (_info.kappa < rhsNorm2)
          _info.stop_tCG = TCGStopReason::reachedTargetResidualKappaLinear;
        else
          _info.stop_tCG = TCGStopReason::reachedTargetResidualThetaSuperLinear;
        break;
      }
      if (residualNorm2 < threshold)
        break;

      z = precond.solve(residual); // approximately solve for "A z = residual"

      RealScalar absOld = absNew;
      absNew            = numext::real(residual.dot(z)); // update the absolute value of r
      RealScalar beta   = absNew / absOld; // calculate the Gram-Schmidt value used to create the new search direction

      e_Pd = beta * (e_Pd + alpha * d_Pd);
      d_Pd = absNew + beta * beta * d_Pd;

      p = z + beta * p; // update search direction
      i++;
    }
    tol_error          = (residualNorm2 / rhsNorm2);
    iters              = i;
    _info.numInnerIter = i;
  }

} // namespace internal

template <typename MatrixType, int UpLo = Lower,
          typename Preconditioner = DiagonalPreconditioner<typename MatrixType::Scalar> >
class TruncatedConjugateGradient;

namespace internal {

  template <typename MatrixType_, int UpLo, typename Preconditioner_>
  struct traits<TruncatedConjugateGradient<MatrixType_, UpLo, Preconditioner_> >
  {
    typedef MatrixType_ MatrixType;
    typedef Preconditioner_ Preconditioner;
  };

} // namespace internal

/**
 * @class TruncatedConjugateGradient
 * \brief Iterative solver for solving linear systems using the truncated conjugate gradient method.
 * \tparam M Type of the matrix A.
 * \tparam upLo Type of the triangular part of the matrix (Lower or Upper or both).
 * \tparam PC Type of the preconditioner.
 */
template <typename M, int upLo, typename PC>
class TruncatedConjugateGradient : public IterativeSolverBase<TruncatedConjugateGradient<M, upLo, PC> >
{
public:
  typedef IterativeSolverBase<TruncatedConjugateGradient> Base;
  TruncatedConjugateGradient(TruncatedConjugateGradient&& other) noexcept
      : Base(std::move(other)),
        algInfo_{other.algInfo_} {}

private:
  using Base::m_error;
  using Base::m_info;
  using Base::m_isInitialized;
  using Base::m_iterations;
  using Base::matrix;
  mutable TCGInfo<typename M::RealScalar> algInfo_;

public:
  using MatrixType     = M;
  using Scalar         = typename MatrixType::Scalar;
  using RealScalar     = typename MatrixType::RealScalar;
  using Preconditioner = PC;

  enum
  {
    UpLo = upLo
  };

public:
  /**
   * \brief Get information about the truncated conjugate gradient algorithm.
   * \return Information about the algorithm.
   */
  TCGInfo<typename MatrixType::RealScalar> getInfo() { return algInfo_; }

  /**
   * \brief Set information about the truncated conjugate gradient algorithm.
   * \param _alginfo Information about the algorithm.
   */
  void setInfo(TCGInfo<typename MatrixType::RealScalar> alginfo) { this->algInfo_ = alginfo; }
  /** Default constructor. */
  TruncatedConjugateGradient()
      : Base() {}

  /** Initialize the solver with matrix \a A for further \c Ax=b solving.
   *
   * This constructor is a shortcut for the default constructor followed
   * by a call to compute().
   *
   * \warning this class stores a reference to the matrix A as well as some
   * precomputed values that depend on it. Therefore, if \a A is changed
   * this class becomes invalid. Call compute() to update it with the new
   * matrix A, or modify a copy of A.
   */
  template <typename MatrixDerived>
  explicit TruncatedConjugateGradient(const EigenBase<MatrixDerived>& A)
      : Base(A.derived()) {}

  ~TruncatedConjugateGradient() {}

  /** \internal */
  template <typename Rhs, typename Dest>
  void _solve_vector_with_guess_impl(const Rhs& b, Dest& x) const {
    typedef typename Base::MatrixWrapper MatrixWrapper;
    typedef typename Base::ActualMatrixType ActualMatrixType;
    enum
    {
      TransposeInput = (!MatrixWrapper::MatrixFree) && (UpLo == (Lower | Upper)) && (!MatrixType::IsRowMajor) &&
                       (!NumTraits<Scalar>::IsComplex)
    };
    typedef std::conditional_t<TransposeInput, Transpose<const ActualMatrixType>, const ActualMatrixType&>
        RowMajorWrapper;
    EIGEN_STATIC_ASSERT(internal::check_implication(MatrixWrapper::MatrixFree, UpLo == (Lower | Upper)),
                        MATRIX_FREE_CONJUGATE_GRADIENT_IS_COMPATIBLE_WITH_UPPER_UNION_LOWER_MODE_ONLY);
    typedef std::conditional_t<UpLo == (Lower | Upper), RowMajorWrapper,
                               typename MatrixWrapper::template ConstSelfAdjointViewReturnType<UpLo>::Type>
        SelfAdjointWrapper;

    m_iterations = Base::maxIterations();
    m_error      = Base::m_tolerance;

    RowMajorWrapper row_mat(matrix());
    internal::truncated_conjugate_gradient(SelfAdjointWrapper(row_mat), b, x, Base::m_preconditioner, m_iterations,
                                           m_error, algInfo_);
    m_info = m_error <= Base::m_tolerance ? Success : NoConvergence;
  }
};

} // end namespace Eigen
