// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file modalanalysishelper.hh
 * \brief Helper for modal analysis related things
 */

#pragma once
#include <dune/vtk/pvdwriter.hh>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/io/vtkwriter.hh>
#include <ikarus/solver/eigenvaluesolver/generalizedeigensolver.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/makeenum.hh>
#include <ikarus/utils/modalanalysis/lumpingschemes.hh>

namespace Ikarus::Dynamics {

/**
 * \brief Factory function to create a lumped mass Assembler from an underlying assembler. You can provide a
 * lumpingscheme similar to those found at \file ikarus/utils/modalanalysis/lumpingschemes.hh. Defaults to row-sum
 * lumping.
 *
 * \tparam AS type of the underlying assembler
 * \tparam LumpingScheme type of the lumping scheme, defaults to row-sum lumping
 * \param assembler the underlying assembler
 * \return auto Assembler Manipulator for mass lumping.
 */
template <Concepts::FlatAssembler AS, typename LumpingScheme = LumpingSchemes::RowSumLumping>
auto makeLumpedFlatAssembler(const std::shared_ptr<AS>& assembler) {
  auto lumpedAssembler = makeAssemblerManipulator(*assembler);
  lumpedAssembler->bind(LumpingScheme{});
  return lumpedAssembler;
}

/**
 * \brief Writes the first nev_ províded eigenmodes from an eigenvalue solver to a vtk file. Each mode is hereby a
 * separate vector quantity.
 *
 * \tparam Eigensolver type of the eigenvalue solver
 * \tparam Assembler type of the assembler
 * \param solver the eigenvalue solver, which provides eigenvectors
 * \param assembler the assembler, which provides the gridview and the finite elements.
 * \param filename Name of the output vtk file
 * \param nev_  optionally specify how many eigenmodes should be written out, defaults to all.
 */
template <Concepts::EigenValueSolver Eigensolver, Concepts::FlatAssembler Assembler>
void writeEigenmodesToVTK(const Eigensolver& solver, std::shared_ptr<Assembler> assembler, const std::string& filename,
                          std::optional<Eigen::Index> nev_ = std::nullopt) {
  Impl::assertNev(nev_, solver.nev());
  auto nev          = nev_.value_or(solver.nev());
  auto eigenvectors = solver.eigenvectors();
  auto basis        = assembler->basis();

  auto writer = Ikarus::Vtk::Writer(assembler);
  for (auto i : Dune::range(nev)) {
    auto evG = assembler->createFullVector(eigenvectors.col(i));
    writer.addInterpolation(std::move(evG), basis, "Eigenmode " + std::to_string(i));
  }
  writer.write(filename);
}

/**
 * \brief
 *
 * \tparam Derived Writes the first nev_ províded eigenmodesin a column-matrix to a paraview collection file
 * (*.pvd). Each mode is hereby a separate vector quantity.
 * \tparam Assembler type of the assembler
 * \param eigenvectors the eigenvectors in a matrix, where the columns represent the eigenmodes.
 * \param assembler the assembler, which provides the gridview and the finite elements.
 * \param filename the filename of the vtk file
 */
template <typename Derived, Concepts::FlatAssembler Assembler>
void writeEigenmodesAsTimeSeries(const Eigen::MatrixBase<Derived>& eigenvectors, std::shared_ptr<Assembler> assembler,
                                 const std::string& filename) {
  auto nev   = eigenvectors.cols();
  auto basis = assembler->basis();

  auto writer    = std::make_shared<decltype(Ikarus::Vtk::Writer(assembler))>(Ikarus::Vtk::Writer(assembler));
  auto pvdWriter = Dune::Vtk::PvdWriter(writer);

  Eigen::VectorXd evG(assembler->size());
  writer->addInterpolation(evG, basis, "Eigenmode");

  for (auto i : Dune::range(nev)) {
    evG = assembler->createFullVector(eigenvectors.col(i));
    pvdWriter.writeTimestep(i + 1, filename, "data", false);
  }
  pvdWriter.write(filename);
}

/**
 * \brief Writes the first nev_ províded eigenmodes from an eigenvalue solver to a paraview collection file
 * (*.pvd). Each mode is hereby written into a separate vtk file and the eigenmodes can then be seen as a timeseries
 * through the pvd file, which aggravates the separate vtk files.
 *
 * \tparam Eigensolver type of the eigenvalue solver
 * \tparam Assembler type of the assembler
 * \param solver the eigenvalue solver, which provides eigenvectors
 * \param assembler the assembler, which provides the gridview and the finite elements.
 * \param filename the filename of the vtk file
 * \param nev_  optionally specify how many eigenmodes should be written out, defaults to all.
 */
template <Concepts::EigenValueSolver Eigensolver, Concepts::FlatAssembler Assembler>
void writeEigenmodesAsTimeSeries(const Eigensolver& solver, std::shared_ptr<Assembler> assembler,
                                 const std::string& filename, std::optional<Eigen::Index> nev_ = std::nullopt) {
  Impl::assertNev(nev_, solver.nev());
  auto nev          = nev_.value_or(solver.nev());
  auto eigenvectors = solver.eigenvectors().leftCols(nev).eval();
  writeEigenmodesAsTimeSeries(eigenvectors, assembler, filename);
}

} // namespace Ikarus::Dynamics
