// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file modalanalysis.hh
 * \brief Implementation of modal analysis
 */

#pragma once
#include <matplot/matplot.h>

#include <dune/vtk/pvdwriter.hh>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/io/vtkwriter.hh>
#include <ikarus/solver/eigenvaluesolver/generalizedeigensolver.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/modalanalysis/lumpingschemes.hh>
#include <ikarus/utils/modalanalysis/modalanalysishelper.hh>

namespace Ikarus::Dynamics {

/**
 * \brief A strongly typed enum class representing the type of result of a modal analysis
 */
MAKE_ENUM(ModalAnalysisResultType, squaredAngularFrequency, angularFrequency, naturalFrequency);

/**
 * \brief Opinionated wrapper class for GeneralizedSymEigenSolver suited for modal analysis
 *
 * \tparam FEC the type of container for finite elements
 * \tparam DV the type of the DirichletValues
 */
template <typename FEC, typename DV>
struct ModalAnalysis
{
  using Assembler     = SparseFlatAssembler<FEC&, DV>;
  using MatrixType    = typename Assembler::MatrixType;
  using ScalarType    = typename MatrixType::Scalar;
  using FERequirement = typename Assembler::FERequirement;
  using FEContainer   = FEC;

  using LumpedAssembler =
      AssemblerManipulator<Assembler, Ikarus::Impl::AssemblerInterfaceHelper<ScalarAssembler, ScalarManipulator>,
                           Ikarus::Impl::AssemblerInterfaceHelper<VectorAssembler, VectorManipulator>,
                           Ikarus::Impl::AssemblerInterfaceHelper<MatrixAssembler, MatrixManipulator>>;
  using Solver = GeneralizedSymEigenSolver<EigenValueSolverType::Spectra, MatrixType>;

  /**
   * \brief Construct a new Modal Analysis object
   *
   * \tparam FES deduced type of the container of the finite elements passed to the constructor
   * \param fes the container of the finite elements
   * \param dv the DirichletValues
   */
  template <typename FES>
  ModalAnalysis(FES&& fes, const DV& dv)
      : fes_(std::forward<FES>(fes)),
        d_(dv.basis().size()),
        dRef_(d_),
        stiffAssembler_(std::make_shared<Assembler>(fes_, dv)),
        massAssembler_(std::make_shared<Assembler>(fes_, dv)) {
    if constexpr (std::remove_cvref_t<decltype(fes.front())>::Traits::useEigenRef)
      req_.insertGlobalSolution(dRef_);
    else
      req_.insertGlobalSolution(d_);

    stiffAssembler_->bind(req_, Ikarus::AffordanceCollections::elastoStatics, Ikarus::DBCOption::Reduced);
    massAssembler_->bind(req_, Ikarus::AffordanceCollections::modalAnalysis, Ikarus::DBCOption::Reduced);
    lumpedMassAssembler_ = makeAssemblerManipulator(*massAssembler_);
  }

  /**
   * \brief Binds a lumping scheme for the mass matrix. The Lumping scheme should have an operator () and modify the
   * matrix mat. For details see already existing lumping schemes at \file ikarus/utils/modalanalysis/lumpingschemes.hh.
   * It resets already bound matrix manipulation functions.
   *
   * \tparam LumpingScheme The type of the lumping scheme, for example one found at \file
   * ikarus/utils/modalanalysis/lumpingschemes.hh.
   * \param ls the instantiated LumpingScheme object (pass either by value or by template definition).
   */
  template <typename LumpingScheme>
  requires(
      std::is_invocable_v<LumpingScheme, Assembler, const FERequirement&, MatrixAffordance, DBCOption, MatrixType&>)
  void bindLumpingScheme(LumpingScheme ls = {}) {
    lumpedMassAssembler_->bind(ls);
  }

  /**
   * \brief Unbinds all former bound lumpingscheme.
   * \remark We have no way to keep track of the applied lumping schemes, therefore we can only unbind all previously
   * bound matrix manipulator functions. In this case however this does not have any unwanted side effects, as the
   * assembler is only used inside the class.p
   */
  void unBindLumpingSchemes() { lumpedMassAssembler_->unbindAllMatrixFunctions(); }

  /**
   * \brief Starts the computation of the eigenvalue solver
   *
   * \param tolerance given tolerance for iterative eigenvalue solving (default: 1e-10)
   * \param maxit givenn maximum iterations for eigenvalue solving (default: 1000)
   * \return true solving was successful
   * \return false solving was not successful
   */
  bool compute(Eigen::Index maxit = 1000, ScalarType tolerance = 1e-10) {
    solver_.emplace(stiffAssembler_, lumpedMassAssembler_);
    return solver_->compute(maxit, tolerance);
  }

  /**
   * \brief Returns the angular frequencies as \f$ \omega  = \sqrt{\lambda} \f$, with \f$ \lambda \f$: eigenvalues from
   * the eigenvalue solver
   */
  Eigen::VectorXd angularFrequencies() {
    assertCompute();
    return squaredAngularFrequencies().cwiseSqrt().eval();
  }

  /**
   * \brief Returns the angular frequencies as \f$ f  = \dfrac{\omega}{2\pi} \f$, with \f$ \omega \f$: angular frequency
   */
  Eigen::VectorXd naturalFrequencies() {
    assertCompute();
    return angularFrequencies() / (2 * std::numbers::pi);
  }

  /**
   * \brief Returns the angular frequencies as \f$ \omega^2  = \lambda \f$, with \f$ \lambda \f$: eigenvalues from
   * the eigenvalue solver
   */
  const Eigen::VectorXd& squaredAngularFrequencies() const {
    assertCompute();
    return solver_->eigenvalues();
  }

  /**
   * \brief Returns the eigenmodes (eigenvectors) of the general eigenvalue problem
   *
   * \return auto matrix with the eigevectors as columns
   */
  const Eigen::MatrixXd& eigenmodes() const { return solver_->eigenvectors(); }

  /**
   * \brief Plots the spectrum of the specified result (defaults to angular frequency) using matplot++.
   *
   * \param resultType specified result type of the modal analysis.
   * \param normalizeModeNumber if true normalizes the output spectrum to [0, 1].
   */
  void plotModalSpectrum(ModalAnalysisResultType resultType = ModalAnalysisResultType::angularFrequency,
                         bool normalizeModeNumber           = false) {
    assertCompute();
    using namespace matplot;
    auto freq = frequencies(resultType);

    auto modeNumbersView = std::ranges::iota_view{1l, nev() + 1} | std::ranges::views::transform([&](auto i) {
                             if (normalizeModeNumber)
                               return static_cast<double>(i) / nev();
                             return static_cast<double>(i);
                           });

    std::vector<double> modeNumbers(modeNumbersView.begin(), modeNumbersView.end());

    auto fig = figure(true);
    auto ax  = fig->add_axes();
    ax->plot(modeNumbers, freq)->line_width(1).color("b");
    ax->xlabel("Mode Number");
    ax->ylabel(toString(resultType));
    ax->title("Modal Spectrum");
    ax->grid(true);

    // Show the figure
    matplot::show();
  }

  /**
   * \brief Writes the first nev_ eigenmodes to a paraview collection file (*.pvd).
   *
   * \param filename filename Name of the output pvd file.
   * \param nev_ optionally specify how many eigenmodes should be written out, defaults to all.
   */
  void writeEigenModes(const std::string& filename, std::optional<Eigen::Index> _nev = std::nullopt) const {
    assertCompute();
    Impl::assertNev(_nev, nev());
    writeEigenmodesAsTimeSeries(solver_.value(), stiffAssembler_, filename, _nev);
  }

  /** \brief Returns the number of eigenvalues of the problem */
  auto nev() const { return solver_->nev(); }

  /** \brief Returns a const reference to the assembler of the stiffness matrix */
  auto& stiffnessAssembler() const { return stiffAssembler_; }

  /** \brief Returns a const reference to the assembler of the (potentially lumped) mass matrix */
  auto& massAssembler() const { return lumpedMassAssembler_; }

private:
  FEContainer fes_;
  FERequirement req_{};
  Eigen::VectorXd d_;
  Eigen::Ref<Eigen::VectorXd> dRef_;
  std::shared_ptr<Assembler> stiffAssembler_;
  std::shared_ptr<Assembler> massAssembler_;
  std::shared_ptr<LumpedAssembler> lumpedMassAssembler_;

  std::optional<Solver> solver_{};

  void assertCompute() const {
    if (not solver_)
      DUNE_THROW(Dune::IOError, "Eigenvalues and -vectors not yet computed, please call compute() first");
  }

  /**
   * \brief Returns one of the above results according to a a specified result type
   *
   * \param rt specified result type of the modal analysis
   */
  auto frequencies(ModalAnalysisResultType rt) {
    if (rt == ModalAnalysisResultType::angularFrequency)
      return angularFrequencies();
    if (rt == ModalAnalysisResultType::naturalFrequency)
      return naturalFrequencies();
    if (rt == ModalAnalysisResultType::squaredAngularFrequency)
      return squaredAngularFrequencies();
    DUNE_THROW(Dune::NotImplemented, "Requested result not implemented");
  }
};

#ifndef DOXYGEN
template <typename FEC, typename DV>
ModalAnalysis(FEC&&, const DV&) -> ModalAnalysis<FEC, DV>;
#endif

} // namespace Ikarus::Dynamics
