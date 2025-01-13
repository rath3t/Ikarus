// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file mixin.hh
 * \brief Implementation of the finite element CRTP mixin class
 * \ingroup  finiteelements
 */

#pragma once

#include <dune/common/tuplevector.hh>

#include <ikarus/finiteelements/fetraits.hh>
#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {
/**
 * \brief CRTP mixin class for finite elements with additional skills.
 *
 * This mixin class is designed for finite elements and provides the ability to include additional skills
 * through the CRTP (Curiously Recurring Template Pattern).
 *
 * @tparam PreFE The base finite element class.
 * @tparam Skills A template parameter pack for additional skills to be mixed into the finite element.
 */
template <typename PreFE, template <typename, typename> class... Skills>
struct FEMixin : public IObserver<NonLinearSolverMessages>, Skills<PreFE, typename PreFE::template FE<Skills...>>...
{
  /**
   * \brief Constructor for the FEMixin class.
   *
   * @param skillsArgs Variadic arguments for initializing the additional skills.
   */
  explicit FEMixin(typename Skills<PreFE, typename PreFE::template FE<Skills...>>::Pre&&... skillsArgs)
      : Skills<PreFE, typename PreFE::template FE<Skills...>>(
            std::forward<typename Skills<PreFE, typename PreFE::template FE<Skills...>>::Pre>(skillsArgs))... {}

  /**
   * \brief Checks if the mixin class has a specific skill.
   *
   * @tparam Skill The skill to check for.
   * @return True if the skill is present, false otherwise.
   */
  template <template <typename, typename> class Skill>
  static consteval bool hasSkill() {
    return Ikarus::traits::hasType<Skill<PreFE, typename PreFE::template FE<Skills...>>,
                                   std::tuple<Skills<PreFE, typename PreFE::template FE<Skills...>>...>>::value;
  }

private:
  template <typename T>
  consteval static auto computeSupportedResultTypes() {
    if constexpr (requires { typename T::SupportedResultTypes; })
      return typename T::SupportedResultTypes();
    else
      return std::tuple();
  }

public:
  /**
   * \brief Type alias for the supported result types by the mixin.
   */
  using SupportedResultTypes =
      decltype(std::tuple_cat(computeSupportedResultTypes<Skills<PreFE, typename PreFE::template FE<Skills...>>>()...));

  template <bool, typename = void>
  struct RequirementType;

  template <typename T>
  struct RequirementType<false, T>
  {
    using type = FERequirements<FESolutions::noSolution, FEParameter::noParameter>;
  };

  template <typename T>
  struct RequirementType<true, T>
  {
    using type = std::common_type_t<typename Skills<PreFE, typename PreFE::template FE<Skills...>>::Requirement...>;
  };

private:
  static constexpr bool requirementDetected =
      Dune::Std::is_detected_v<std::common_type_t,
                               typename Skills<PreFE, typename PreFE::template FE<Skills...>>::Requirement...>;
  static_assert(requirementDetected or sizeof...(Skills) == 0, "The skills must have a common fe requirement type.");

public:
  using Traits                  = PreFE::Traits;
  using Requirement             = RequirementType<requirementDetected>::type;
  using LocalView               = typename Traits::LocalView;
  static constexpr int worldDim = Traits::worlddim;

  template <StrainTags ES>
  static constexpr bool hasEAS = hasSkill<EnhancedAssumedStrainsPre<ES>::template Skill>();

  /**
   * \brief Create a Requirement object.
   *
   * @return The created Requirement object.
   */
  static auto createRequirement() { return Requirement(); }

  /**
   * \brief Calculate the scalar value associated with the given Requirement.
   *
   * \tparam ScalarType The scalar type for the calculation.
   * \param req The Requirement object specifying the requirements for the calculation.
   * \return The calculated scalar value.
   */
  friend auto calculateScalar(const FEMixin& self, const Requirement& req, ScalarAffordance affordance) {
    return self.template calculateScalarImpl<double>(req, affordance);
  }

  /**
   * \brief Calculate the vector associated with the given Requirement.
   *
   * \tparam ScalarType The scalar type for the calculation.
   * \param req The Requirement object specifying the requirements for the calculation.
   * \param force The vector to store the calculated result.
   */
  friend void calculateVector(const FEMixin& self, const Requirement& req, VectorAffordance affordance,
                              typename Traits::template VectorType<> force) {
    self.template calculateVectorImpl<double>(req, affordance, force);
  }

  /**
   * \brief Calculate the matrix associated with the given Requirement.
   *
   * \tparam ScalarType The scalar type for the calculation.
   * \param req The Requirement object specifying the requirements for the calculation.
   * \param K The matrix to store the calculated result.
   */
  friend void calculateMatrix(const FEMixin& self, const Requirement& req, MatrixAffordance affordance,
                              typename Traits::template MatrixType<> K) {
    self.template calculateMatrixImpl<double>(req, affordance, K);
  }

  using Skills<PreFE, typename PreFE::template FE<Skills...>>::calculateAtImpl...;

  /**
   * \brief Calculate the element values at a specific location for a given ResultType.
   *
   * @tparam RT The ResultType to calculate.
   * @param req The Requirement object specifying the requirements for the calculation.
   * @param local The local coordinates where the calculation is performed.
   * @return The calculated result as specified by the ResultType.
   */
  template <template <typename, int, int> class RT>
  requires requires(FEMixin m, const Requirement& req, const Dune::FieldVector<double, Traits::mydim>& local) {
    m.template calculateAtImpl<RT>(req, local, Dune::PriorityTag<10>());
  }
  auto calculateAt(const Requirement& req, const Dune::FieldVector<double, Traits::mydim>& local) const {
    return this->template calculateAtImpl<RT>(req, local, Dune::PriorityTag<10>());
  }

private:
  template <typename Sk>
  auto invokeBind() {
    if constexpr (requires { this->Sk::bindImpl(); })
      Sk::bindImpl();
  }

  static constexpr bool implementsCalculateScalarImpl =
      (requires(FEMixin m, const Requirement& par, ScalarAffordance affordance,
                const std::optional<std::reference_wrapper<const Eigen::VectorX<double>>>& dx) {
        m.Skills<PreFE, typename PreFE::template FE<Skills...>>::calculateScalarImpl(par, affordance, dx);
      } ||
       ...);

public:
  /**
   * \brief Call all bind functions if the skill implements it
   */
  void bind() { (invokeBind<Skills<PreFE, typename PreFE::template FE<Skills...>>>(), ...); }

  /**
   * \brief Calculate the scalar value in each skill and joins them by `+`.
   *
   * @tparam ScalarType The scalar type for the calculation.
   * @param par The Requirement object specifying the requirements for the calculation.
   * @param dx Optional vector used in the calculation.
   * @return The calculated scalar value.
   */
  template <typename ScalarType = double>
  requires implementsCalculateScalarImpl
  auto calculateScalarImpl(
      const Requirement& par, ScalarAffordance affordance,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    return (Skills<PreFE, typename PreFE::template FE<Skills...>>::template calculateScalarImpl<ScalarType>(
                par, affordance, dx) +
            ... + ScalarType{0});
  }

private:
  static constexpr bool implementsCalculateVectorImpl =
      (requires(FEMixin m, const Requirement& par, VectorAffordance affordance,
                typename Traits::template VectorType<double> force,
                const std::optional<std::reference_wrapper<const Eigen::VectorX<double>>>& dx) {
        m.Skills<PreFE, typename PreFE::template FE<Skills...>>::calculateVectorImpl(par, affordance, force, dx);
      } ||
       ...);

public:
  /**
   * \brief Calculate the vector for each skill
   *
   * @tparam ScalarType The scalar type for the calculation.
   * @param par The Requirement object specifying the requirements for the calculation.
   * @param force The vector to store the calculated result.
   * @param dx Optional vector used in the calculation.
   */
  template <typename ScalarType>
  requires implementsCalculateVectorImpl
  void calculateVectorImpl(
      const Requirement& par, VectorAffordance affordance, typename Traits::template VectorType<ScalarType> force,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    (Skills<PreFE, typename PreFE::template FE<Skills...>>::template calculateVectorImpl<ScalarType>(par, affordance,
                                                                                                     force, dx),
     ...);
  }

private:
  static constexpr bool implementsCalculateMatrixImpl =
      (requires(FEMixin m, const Requirement& par, MatrixAffordance affordance,
                typename Traits::template MatrixType<double> K,
                const std::optional<std::reference_wrapper<const Eigen::VectorX<double>>>& dx) {
        m.Skills<PreFE, typename PreFE::template FE<Skills...>>::calculateMatrixImpl(par, affordance, K, dx);
      } ||
       ...);

public:
  /**
   * \brief Calculate the matrix for each skill
   *
   * @tparam ScalarType The scalar type for the calculation.
   * @param par The Requirement object specifying the requirements for the calculation.
   * @param K The matrix to store the calculated result.
   * @param dx Optional vector used in the calculation.
   */
  template <typename ScalarType>
  requires implementsCalculateMatrixImpl
  void calculateMatrixImpl(
      const Requirement& par, MatrixAffordance affordance, typename Traits::template MatrixType<ScalarType> K,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    (Skills<PreFE, typename PreFE::template FE<Skills...>>::template calculateMatrixImpl<ScalarType>(par, affordance, K,
                                                                                                     dx),
     ...);
  }

private:
  template <typename Sk>
  auto invokeUpdateState(NonLinearSolverMessages message, const Requirement& par,
                         const std::remove_reference_t<typename Traits::template VectorType<>>& correction) const {
    if constexpr (requires { Sk::updateStateImpl(message, par, correction); })
      Sk::updateStateImpl(message, par, correction);
  }

public:
  void updateImpl(NonLinearSolverMessages message, Eigen::VectorXd& vec, const Eigen::VectorXd& dx) {
    auto req = Requirement();
    req.insertGlobalSolution(vec);
    updateState(message, req, dx);
  }

  /**
   * \brief  Call all updateStateImpl functions if the skill implements it.
   *
   * \details Update the state variables related to a particular skill.
   *
   * \param req The Requirement object specifying the requirements for the update itself.
   * \param correction A correction vector (for example, the displacement increment) based on which the state variables
   * are to be updated.
   */
  void updateState(NonLinearSolverMessages message, const Requirement& par,
                   const std::remove_reference_t<typename Traits::template VectorType<>>& correction) const {
    (invokeUpdateState<Skills<PreFE, typename PreFE::template FE<Skills...>>>(message, par, correction), ...);
  }

protected:
  /**
   * \brief Get a reference to the underlying finite element object.
   *
   * @return A reference to the underlying finite element object.
   */
  const auto& underlying() const { return static_cast<const typename PreFE::template FE<Skills...>&>(*this); }
  /**
   * \brief Get a reference to the underlying finite element object.
   *
   * @return A reference to the underlying finite element object.
   */
  auto& underlying() { return static_cast<typename PreFE::template FE<Skills...>&>(*this); }
};

/**
 * \brief Struct representing a collection of skills.
 *
 * @tparam ARGS Variadic template parameters representing the skills.
 */
template <typename... ARGS>
struct Skills
{
  using Args = std::tuple<ARGS...>;
  Args args;
};

/**
 * \brief Function to create a Skills instance with the given skills.
 *
 * @tparam Args Variadic template parameters representing the skills.
 * @param args Variadic arguments representing the skills.
 * @return A Skills instance containing the specified skills.
 */
template <typename... Args>
auto skills(const Args&... args) {
  return Skills<std::remove_cvref_t<Args>...>{std::forward_as_tuple(std::remove_cvref_t<Args>(args)...)};
}

/**
 * \brief Function to merge two Skills instances.
 *
 * @tparam Args1 Variadic template parameters representing the skills of the first instance.
 * @tparam Args2 Variadic template parameters representing the skills of the second instance.
 * @param sk1 The first Skills instance.
 * @param sk2 The second Skills instance.
 * @return A new Skills instance containing the merged skills.
 */
template <typename... Args1, typename... Args2>
auto merge(const Skills<Args1...>& sk1, const Skills<Args2...>& sk2) {
  return Skills<std::remove_cvref_t<Args1>..., std::remove_cvref_t<Args2>...>{std::tuple_cat(sk1.args, sk2.args)};
}

} // namespace Ikarus
