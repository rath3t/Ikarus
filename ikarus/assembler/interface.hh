// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file interface.hh
 * \brief Defines the interface for flat, scalar, vector and matrix assemblers
 */

#pragma once
#include <ranges>

#include <dune/common/referencehelper.hh>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/finiteelements/fehelper.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/dirichletvalues.hh>

namespace Ikarus {
/**
 * \class FlatAssemblerBase
 * \brief The FlatAssemblerBase takes care of common subtasks done by flat assemblers.
 * \ingroup assembler
 *
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 */
template <typename FEC, typename DV>
class FlatAssemblerBase
{
public:
  using FEContainerRaw = std::remove_cvref_t<FEC>; ///< Type of the raw finite element container.
  using FERequirement  = typename FEContainerRaw::value_type::Requirement;
  ///< Type of the finite element requirement.
  using GlobalIndex = typename FEContainerRaw::value_type::GlobalIndex; ///< Type of the global index.
  using Basis       = typename DV::Basis;                               ///< Type of the basis.
  using GridView    = typename Basis::GridView;                         ///< Type of the grid view.
  using FEContainer = FEC;                                              ///< Type of the finite element container.
  ///< Type of the finite element container (reference or by value).
  using DirichletValuesType = DV;                          ///< Type of the Dirichlet values.
  using SizeType = typename DirichletValuesType::SizeType; ///< size_type of the container storing Dirichlet flags
  using AffordanceCollectionType = AffordanceCollection<ScalarAffordance, VectorAffordance, MatrixAffordance>;

  /**
   * \brief Constructor for FlatAssemblerBase.
   *
   * \param fes Finite element container.
   * \param dirichletValues Reference to Dirichlet values.
   */
  template <typename FEContainer_ = FEContainer, typename DirichletValuesType_ = DirichletValuesType>
  FlatAssemblerBase(FEContainer_&& fes, DirichletValuesType_&& dirichletValues)
      : feContainer_{std::forward<FEContainer_>(fes)},
        dirichletValues_{std::forward<DirichletValuesType_>(dirichletValues)} {
    constraintsBelow_.reserve(dirichletValues_.size());
    size_t counter = 0;
    for (auto iv : std::ranges::iota_view{decltype(dirichletValues_.size())(0), dirichletValues_.size()}) {
      constraintsBelow_.emplace_back(counter);
      if (dirichletValues_.isConstrained(iv))
        ++counter;
    }
    fixedDofs_ = dirichletValues_.fixedDOFsize();
  }

  /**
   * \brief Returns the size of the free degrees of freedom, which are not fixed by a Dirichlet boundary condition.
   * \return Size of the reduced degrees of freedom.
   */
  size_t reducedSize() { return dirichletValues_.size() - fixedDofs_; }

  /**
   * \brief Returns the size of nodes, i.e., the number of degrees of freedom.
   * \return Size of the degrees of freedom.
   */
  size_t size() { return dirichletValues_.size(); }

  /**
   * \brief Creates the full-sized vector of size Dof and inserts the values of a reduced vector at the "free"
   * degrees of freedom and writes a zero for the fixed doffs.
   *
   * \param reducedVector Reference to the reduced vector.
   * \return Eigen::VectorXd The full-sized vector.
   */
  Eigen::VectorXd createFullVector(Eigen::Ref<const Eigen::VectorXd> reducedVector);

  /**
   * \brief Returns the container of finite elements.
   * \return Reference to the finite element container.
   */
  auto& finiteElements() const { return Dune::resolveRef(feContainer_); }

  /**
   * \brief Returns the dirichlet value object.
   * \return Reference to the dirichlet value object.
   */
  const auto& dirichletValues() const { return Dune::resolveRef(dirichletValues_); }

  /**
   * \brief Returns the basis object
   * \return Reference to the basis
   */
  const auto& basis() const { return Dune::resolveRef(dirichletValues_.basis()); }

  /**
   * \brief Returns the gridView object.
   * \return Reference to the gridView object.
   */
  const auto& gridView() const { return Dune::resolveRef(dirichletValues_.basis().gridView()); }

  /**
   * \brief Returns the number of constraints below a given degrees of freedom index.
   *
   * \param i Index of the degree of freedom.
   * \return Number of constraints below the given index.
   */
  [[nodiscard]] size_t constraintsBelow(SizeType i) const { return constraintsBelow_[i]; }

  /**
   * \brief Returns true if a given degree of freedom is fixed by a Dirichlet boundary condition.
   *
   * \param i Index of the degree of freedom.
   * \return True if the degree of freedom is fixed; false otherwise.
   */
  [[nodiscard]] bool isConstrained(SizeType i) const { return dirichletValues_.isConstrained(i); }

  /**
   * \brief Coarse estimate of node connectivity, i.e., this relates to the bandwidth of a sparse matrix.
   * This estimate overestimates the real connectivity and should only be used for allocating vectors.
   *
   * \return Size_t Coarse estimate of node connectivity.
   */
  [[nodiscard]] size_t estimateOfConnectivity() const {
    return dirichletValues_.basis().gridView().size(GridView::dimension) * 8;
  }

  /**
   * \brief Binds the assembler to a set of finite element requirement and affordance.
   *
   * \param req Reference to the finite element requirement.
   * \param affordanceCollection The affordance collection
   */
  void bind(const FERequirement& req, AffordanceCollectionType affordanceCollection,
            DBCOption dbcOption = DBCOption::Full) {
    req_         = std::make_optional<FERequirement>(req);
    affordances_ = std::make_optional<AffordanceCollectionType>(affordanceCollection);
    dBCOption_   = std::make_optional<DBCOption>(dbcOption);
  }

  /**
   * \brief Binds the assembler to a  finite element requirement.
   *
   * \param req Reference to the finite element requirement.
   */
  void bind(const FERequirement& req) { req_ = std::make_optional<FERequirement>(req); }

  /**
   * \brief Binds the assembler to an affordance collection.
   *
   * \param affordanceCollection The affordance collection
   */
  void bind(AffordanceCollectionType affordanceCollection) {
    affordances_ = std::make_optional<AffordanceCollectionType>(affordanceCollection);
  }

  /**
   * \brief Binds the assembler to an affordance collection.
   *
   * \param dbcOption The EnforcingDBC option
   */
  void bind(DBCOption dbcOption) { dBCOption_ = std::make_optional<DBCOption>(dbcOption); }

  /**
   * \brief Returns true if the assembler is bound to a finite element requirement and affordance.
   *
   * \return True if the assembler is bound; false otherwise.
   */
  [[nodiscard]]
  bool bound() const {
    if (boundToRequirement() and boundToAffordanceCollection() and boundToDBCOption())
      return true;
    else
      DUNE_THROW(Dune::InvalidStateException, "The assembler is not bound to a requirement, affordance or dBCOption.");
  }

  /**
   * \brief Returns true if the assembler is bound to a finite element requirement.
   *
   * \return True if the assembler is bound; false otherwise.
   */
  [[nodiscard]]
  bool boundToRequirement() const {
    return req_.has_value();
  }

  /**
   * \brief Returns true if the assembler is bound to an affordance collection.
   *
   * \return True if the assembler is bound; false otherwise.
   */
  [[nodiscard]]
  bool boundToAffordanceCollection() const {
    return affordances_.has_value();
  }

  /**
   * \brief Returns true if the assembler is bound to an affordance collection.
   *
   * \return True if the assembler is bound; false otherwise.
   */
  [[nodiscard]]
  bool boundToDBCOption() const {
    return dBCOption_.has_value();
  }

  /**
   * \brief Returns the requirement.
   *
   */
  FERequirement& requirement() {
    if (req_.has_value())
      return req_.value();
    else
      DUNE_THROW(Dune::InvalidStateException, "The requirement can only be obtained after binding");
  }

  /**
   * \brief Returns the affordance.
   *
   */
  AffordanceCollectionType affordanceCollection() const {
    if (affordances_.has_value())
      return affordances_.value();
    else
      DUNE_THROW(Dune::InvalidStateException, "The affordance can only be obtained after binding");
  }

  /**
   * \brief Returns the dirichlet boundary condition enforcement option.
   *
   */
  DBCOption dBCOption() const {
    if (dBCOption_.has_value())
      return dBCOption_.value();
    else
      DUNE_THROW(Dune::InvalidStateException, "The dBCOption can only be obtained after binding");
  }

private:
  FEContainer feContainer_;
  DirichletValuesType dirichletValues_;
  std::optional<FERequirement> req_;
  std::optional<AffordanceCollectionType> affordances_;
  std::vector<size_t> constraintsBelow_{};
  size_t fixedDofs_{};
  std::optional<DBCOption> dBCOption_;
};

#ifndef DOXYGEN
/// If we have an rvalue reference then we want to store the type on our own. Otherwise, we just store references
template <class FEV, class DirichletValuesType>
FlatAssemblerBase(const FEV& fes,
                  const DirichletValuesType& dirichletValues) -> FlatAssemblerBase<FEV, DirichletValuesType>;
#endif

/**
 * \brief The ScalarAssembler provides an interface for an assembler that assembles scalar quantities.
 * \tparam SA Type of the scalar assembler
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 * \tparam ST Type of the scalar quantity
 */
template <typename SA, typename FEC, typename DV, typename ST>
class ScalarAssembler
{
public:
  using ScalarAssemblerType = SA;
  using ScalarType          = ST;
  using DirichletValuesType = DV;
  using FEContainer         = FEC;
  using FEContainerRaw      = std::remove_cvref_t<FEC>; ///< Type of the raw finite element container.
  using FERequirement       = typename FEContainerRaw::value_type::Requirement;

  /**
   * \brief Calculates the scalar quantity requested by feRequirements and affordance.
   *
   * \param feRequirements Reference to the finite element requirements.
   * \param affordance The scalar affordance
   * \return Const reference to the calculated scalar quantity.
   */
  const ScalarType& scalar(const FERequirement& feRequirements, ScalarAffordance affordance) {
    return underlying().getScalarImpl(feRequirements, affordance);
  }

  /**
   * \brief Calculates the scalar quantity requested by the bound feRequirements and returns a reference.
   *
   * \return Const reference to the calculated scalar quantity.
   */
  const ScalarType& scalar() {
    return underlying().getScalarImpl(underlying().requirement(),
                                      underlying().affordanceCollection().scalarAffordance());
  }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const ScalarAssemblerType&>(*this); }
  auto& underlying() { return static_cast<ScalarAssemblerType&>(*this); }
};

/**
 * \brief The VectorAssembler provides an interface for an assembler that assembles vector quantities.
 * \tparam VA Type of the vector assembler
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 * \tparam VT Type of the vector quantity
 */
template <typename VA, typename FEC, typename DV, typename VT>
class VectorAssembler
{
public:
  using VectorAssemblerType = VA;
  using VectorType          = VT;

  using FEContainerRaw = std::remove_cvref_t<FEC>; ///< Type of the raw finite element container.
  using FERequirement  = typename FEContainerRaw::value_type::Requirement;
  ///< Type of the finite element requirement.
  using GlobalIndex = typename FEContainerRaw::value_type::GlobalIndex; ///< Type of the global index.

  using DirichletValuesType = DV;
  using FEContainer         = FEC;

  /**
     * \brief Calculates the vectorial quantity requested by the  feRequirements and the affordance.
     Depending on the requested DBCOption, the raw, reduced or full vector is returned.
      Raw means the degrees of freedom associated with dirichlet boundary conditions are not changed.
      Full means that degrees of freedom associated with dirichlet boundary conditions are set to zero in the vector.
      Reduced means that degrees of freedom associated with dirichlet boundary conditions are removed and the returned
     vector has reduced size.
     *
     * \param feRequirements Reference to the finite element requirements.
     * \param affordance The vector affordance
     * \param dbcOption The DBCOption
     * \return Const reference to the calculated vectorial quantity.
     */
  const VectorType& vector(const FERequirement& feRequirements, VectorAffordance affordance,
                           DBCOption dbcOption = DBCOption::Full) {
    if (dbcOption == DBCOption::Raw) {
      return underlying().getRawVectorImpl(feRequirements, affordance);
    } else if (dbcOption == DBCOption::Reduced) {
      return underlying().getReducedVectorImpl(feRequirements, affordance);
    } else if (dbcOption == DBCOption::Full) {
      return underlying().getVectorImpl(feRequirements, affordance);
    }
    __builtin_unreachable();
  }

  /**
 * \brief Calculates the vectorial quantity requested by the bound feRequirements and the affordance.
 Depending on the requested DBCOption, the raw, reduced or full vector is returned.
  Raw means the degrees of freedom associated with dirichlet boundary conditions are not changed.
  Full means that degrees of freedom associated with dirichlet boundary conditions are set to zero in the vector.
  Reduced means that degrees of freedom associated with dirichlet boundary conditions are removed and the returned
 vector has reduced size.
 * \param dbcOption The DBCOption
 * \return Const reference to the calculated vectorial quantity.
 */
  const VectorType& vector(DBCOption dbcOption) {
    return vector(underlying().requirement(), underlying().affordanceCollection().vectorAffordance(), dbcOption);
  }

  /**
* \brief Calculates the vectorial quantity requested by the bound feRequirements,  the affordance and the
dBCOption. Depending on the DBCOption, the raw, reduced or full vector is returned. Raw
means the degrees of freedom associated with dirichlet boundary conditions are not changed. Full means that degrees of
freedom associated with dirichlet boundary conditions are set to zero in the vector. Reduced means that degrees of
freedom associated with dirichlet boundary conditions are removed and the returned vector has reduced size.
* \return Const reference to the calculated vectorial quantity.
*/
  const VectorType& vector() { return vector(underlying().dBCOption()); }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const VectorAssemblerType&>(*this); }
  auto& underlying() { return static_cast<VectorAssemblerType&>(*this); }
};

/**
 * \brief The MatrixAssembler provides an interface for an assembler that assembles matrix quantities.
 * \tparam MA Type of the matrix assembler
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 * \tparam MT Type of the matrix quantity
 */
template <typename MA, typename FEC, typename DV, typename MT>
class MatrixAssembler
{
public:
  using MatrixAssemblerType = MA;
  using MatrixType          = MT;

  using FEContainerRaw = std::remove_cvref_t<FEC>; ///< Type of the raw finite element container.
  using FERequirement  = typename FEContainerRaw::value_type::Requirement;
  ///< Type of the finite element requirement.
  using GlobalIndex = typename FEContainerRaw::value_type::GlobalIndex; ///< Type of the global index.

  using DirichletValuesType = DV;
  using FEContainer         = FEC;

  /**
   * \brief Calculates the matrix quantity requested by feRequirements and the affordance.
   * For DBCOption::Full a zero is written on fixed degrees of freedom rows and columns, and a one is written
   * on the diagonal. For DBCOption::Raw the untouched matrix is returned.
   * For DBCOption::Reduced the matrix is reduced in size by removing the fixed degrees of freedom.

    \param feRequirements Reference to the finite element requirements.
   * \param affordance The matrix affordance
   * \param dbcOption The DBCOption
   * \return Const reference to the modified sparse matrix quantity.
   */
  const MatrixType& matrix(const FERequirement& feRequirements, MatrixAffordance affordance,
                           DBCOption dbcOption = DBCOption::Full) {
    if (dbcOption == DBCOption::Raw) {
      return underlying().getRawMatrixImpl(feRequirements, affordance);
    } else if (dbcOption == DBCOption::Reduced) {
      return underlying().getReducedMatrixImpl(feRequirements, affordance);
    } else if (dbcOption == DBCOption::Full) {
      return underlying().getMatrixImpl(feRequirements, affordance);
    }
    __builtin_unreachable();
  }

  /**
   * \brief Calculates the matrix quantity requested by the bound feRequirements and the affordance.
   * \see const Eigen::SparseMatrix<double>& matrix(const FERequirement& feRequirements,MatrixAffordance affordance,
   DBCOption dbcOption)

   * \param dbcOption The DBCOption
   * \return Const reference to the modified sparse matrix quantity.
   */
  const MatrixType& matrix(DBCOption dbcOption) {
    return matrix(underlying().requirement(), underlying().affordanceCollection().matrixAffordance(), dbcOption);
  }

  /**
 * \brief Calculates the matrix quantity requested by the bound feRequirements, the affordance and the
dBCOption.
 * \see const Eigen::SparseMatrix<double>& matrix(const FERequirement& feRequirements,MatrixAffordance affordance,
 DBCOption dbcOption)

 * \return Const reference to the modified sparse matrix quantity.
 */
  const MatrixType& matrix() { return matrix(underlying().dBCOption()); }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const MatrixAssemblerType&>(*this); }
  auto& underlying() { return static_cast<MatrixAssemblerType&>(*this); }
};
} // namespace Ikarus
