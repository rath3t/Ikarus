// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearelastic.hh
 * \brief Definition of the NonLinearElastic class for finite element mechanics computations.
 * \ingroup  mechanics
 */

#pragma once

#if HAVE_DUNE_LOCALFEFUNCTIONS
  #include <dune/fufem/boundarypatch.hh>
  #include <dune/geometry/quadraturerules.hh>
  #include <dune/geometry/type.hh>
  #include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
  #include <dune/localfefunctions/expressions/greenLagrangeStrains.hh>
  #include <dune/localfefunctions/impl/standardLocalFunction.hh>
  #include <dune/localfefunctions/manifolds/realTuple.hh>

  #include <ikarus/finiteelements/febase.hh>
  #include <ikarus/finiteelements/feconfiguration.hh>
  #include <ikarus/finiteelements/fehelper.hh>
  #include <ikarus/finiteelements/mechanics/loads.hh>
  #include <ikarus/finiteelements/mechanics/materials/tags.hh>
  #include <ikarus/finiteelements/physicshelper.hh>
  #include <ikarus/utils/defaultfunctions.hh>
  #include <ikarus/utils/eigendunetransformations.hh>
  #include <ikarus/utils/linearalgebrahelper.hh>

namespace Ikarus {

template <typename PreFE, typename FE, typename PRE>
class NonLinearElastic;

/**
 * \brief A PreFE struct for non-linear elastic elements.
 * \tparam MAT Type of the material.
 */
template <typename MAT>
struct NonLinearElasticPre
{
  using Material = MAT;
  MAT material;

  template <typename PreFE, typename FE>
  using Skill = NonLinearElastic<PreFE, FE, NonLinearElasticPre>;
};

/**
 * \brief NonLinearElastic class represents a non-linear elastic finite element.
 *
 * \ingroup mechanics
 *
 * \tparam PreFE The type of the total pre finite element.
 * \tparam FE The type of the finite element.
 * \tparam PRE The type of the non-linear elastic pre finite element.
 */
template <typename PreFE, typename FE, typename PRE>
class NonLinearElastic : public ResultTypeBase<ResultTypes::PK2Stress, ResultTypes::PK2StressFull>
{
public:
  using Traits    = PreFE::Traits;
  using Basis     = typename Traits::Basis;
  using FlatBasis = typename Traits::FlatBasis;
  using Configuration =
      FEConfigurationFactory<FESolutions::displacement, FEParameter::loadfactor, Traits::useEigenRef>::type;
  using LocalView = typename Traits::LocalView;
  using Geometry  = typename Traits::Geometry;
  using GridView  = typename Traits::GridView;
  using Element   = typename Traits::Element;
  using Material  = PRE::Material;
  using Pre       = PRE;

  using LocalBasisType = decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis());

  static constexpr int myDim       = Traits::mydim;
  static constexpr auto strainType = StrainTags::greenLagrangian;

  /**
   * \brief Constructor for the NonLinearElastic class.
   * \param pre The pre fe
   */
  explicit NonLinearElastic(const Pre& pre)
      : mat_{pre.material} {}

protected:
  /**
   * \brief A helper function to bind the local view to the element.
   */
  void bindImpl() {
    const auto& localView = underlying().localView();
    const auto& element   = localView.element();
    auto& firstChild      = localView.tree().child(0);
    const auto& fe        = firstChild.finiteElement();
    geo_                  = std::make_shared<const Geometry>(element.geometry());
    numberOfNodes_        = fe.size();
    order_                = 2 * (fe.localBasis().order());
    localBasis_           = Dune::CachedLocalBasis(fe.localBasis());
    if constexpr (requires { element.impl().getQuadratureRule(order_); })
      if (element.impl().isTrimmed())
        localBasis_.bind(element.impl().getQuadratureRule(order_), Dune::bindDerivatives(0, 1));
      else
        localBasis_.bind(Dune::QuadratureRules<double, myDim>::rule(element.type(), order_),
                         Dune::bindDerivatives(0, 1));
    else
      localBasis_.bind(Dune::QuadratureRules<double, myDim>::rule(element.type(), order_), Dune::bindDerivatives(0, 1));
  }

public:
  /**
   * \brief Get the displacement function for the given Requirement.
   *
   * \tparam ScalarType The scalar type for the displacement function.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \return A StandardLocalFunction representing the displacement function.
   */
  template <typename ScalarType = double>
  auto displacementFunction(
      const Configuration& par,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    const auto& d = par.globalSolution();
    auto disp     = Ikarus::FEHelper::localSolutionBlockVector<Traits>(d, underlying().localView(), dx);
    Dune::StandardLocalFunction uFunction(localBasis_, disp, geo_);
    return uFunction;
  }

  /**
   * \brief The strain function for the given Requirement.
   *
   * \tparam ScalarType The scalar type for the strain function.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \return The strain function calculated using greenLagrangeStrains.
   */
  template <typename ScalarType = double>
  inline auto strainFunction(
      const Configuration& par,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    return Dune::greenLagrangeStrains(displacementFunction(par, dx));
  }

  /**
   * \brief Get the material tangent for the given strain.
   *
   * \tparam ScalarType The scalar type for the material and strain.
   * \tparam strainDim The dimension of the strain vector.
   * \tparam voigt Flag indicating whether to use Voigt notation.
   * \param strain The strain vector.
   * \return The material tangent calculated using the material's tangentModuli function.
   */
  template <typename ScalarType, int strainDim, bool voigt = true>
  auto materialTangent(const Eigen::Vector<ScalarType, strainDim>& strain) const {
    return material<ScalarType>().template tangentModuli<strainType, voigt>(strain);
  }

  /**
   * \brief Get the internal energy for the given strain.
   *
   * \tparam ScalarType The scalar type for the material and strain.
   * \tparam strainDim The dimension of the strain vector.
   * \param strain The strain vector.
   * \return The internal energy calculated using the material's storedEnergy function.
   */
  template <typename ScalarType, int strainDim>
  auto getInternalEnergy(const Eigen::Vector<ScalarType, strainDim>& strain) const {
    return material<ScalarType>().template storedEnergy<strainType>(strain);
  }

  /**
   * \brief Get the stress for the given strain.
   *
   * \tparam ScalarType The scalar type for the material and strain.
   * \tparam strainDim The dimension of the strain vector.
   * \tparam voigt A boolean indicating whether to use the Voigt notation for stress.
   * \param strain The strain vector.
   * \return The stress vector calculated using the material's stresses function.
   */
  template <typename ScalarType, int strainDim, bool voigt = true>
  auto getStress(const Eigen::Vector<ScalarType, strainDim>& strain) const {
    return material<ScalarType>().template stresses<strainType, voigt>(strain);
  }

  const Geometry& geometry() const { return *geo_; }
  [[nodiscard]] size_t numberOfNodes() const { return numberOfNodes_; }
  [[nodiscard]] int order() const { return order_; }

  template <typename ScalarType = double>
  decltype(auto) material() const {
    if constexpr (Concepts::AutodiffScalar<ScalarType>)
      return mat_.template rebind<ScalarType>();
    else
      return mat_;
  }

public:
  /**
   * \brief Calculates a requested result at a specific local position.
   *
   * \param req The Requirement object holding the global solution.
   * \param local Local position vector.
   * \return calculated result
   *
   * \tparam RT The type representing the requested result.
   */
  template <template <typename, int, int> class RT>
  requires(supportsResultType<RT>())
  auto calculateAtImpl(const Configuration& req, const Dune::FieldVector<double, Traits::mydim>& local,
                       Dune::PriorityTag<1>) const {
    using namespace Dune::DerivativeDirections;

    using RTWrapper = ResultWrapper<RT<typename Traits::ctype, myDim, Traits::worlddim>, ResultShape::Vector>;
    if constexpr (isSameResultType<RT, ResultTypes::PK2Stress> or isSameResultType<RT, ResultTypes::PK2StressFull>) {
      const auto uFunction = displacementFunction(req);
      const auto H         = uFunction.evaluateDerivative(local, Dune::wrt(spatialAll), Dune::on(gridElement));
      const auto E         = (0.5 * (H.transpose() + H + H.transpose() * H)).eval();

      decltype(auto) mat = [&]() {
        if constexpr (isSameResultType<RT, ResultTypes::PK2StressFull> and requires { mat_.underlying(); })
          return mat_.underlying();
        else
          return mat_;
      }();
      return RTWrapper{mat.template stresses<StrainTags::greenLagrangian>(enlargeIfReduced<Material>(toVoigt(E)))};
    }
  }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const FE&>(*this); }
  auto& underlying() { return static_cast<FE&>(*this); }
  std::shared_ptr<const Geometry> geo_;
  Dune::CachedLocalBasis<std::remove_cvref_t<LocalBasisType>> localBasis_;
  Material mat_;
  size_t numberOfNodes_{0};
  int order_{};

protected:
  /**
   * \brief Calculate the matrix associated with the given Requirement.
   *
   * \tparam ScalarType The scalar type for the calculation.
   * \param par The Requirement object specifying the requirements for the calculation.
   * \param K The matrix to store the calculated result.
   */
  template <typename ScalarType>
  void calculateMatrixImpl(
      const Configuration& par, const MatrixAffordance& affordance, typename Traits::template MatrixType<> K,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto uFunction = displacementFunction(par, dx);
    const auto eps       = strainFunction(par, dx);
    for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
      const double intElement = geo_->integrationElement(gp.position()) * gp.weight();
      const auto EVoigt       = (eps.evaluate(gpIndex, on(gridElement))).eval();
      const auto u            = (uFunction.evaluate(gpIndex, on(gridElement))).eval();
      const auto C            = materialTangent(EVoigt);

      const auto stresses = getStress(EVoigt);
      for (size_t i = 0; i < numberOfNodes_; ++i) {
        const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
        for (size_t j = 0; j < numberOfNodes_; ++j) {
          const auto bopJ = eps.evaluateDerivative(gpIndex, wrt(coeff(j)), on(gridElement));
          const auto kgIJ = eps.evaluateDerivative(gpIndex, wrt(coeff(i, j)), along(stresses), on(gridElement));
          K.template block<myDim, myDim>(i * myDim, j * myDim) += (bopI.transpose() * C * bopJ + kgIJ) * intElement;
        }
      }
    }
  }

  template <typename ScalarType>
  auto calculateScalarImpl(const Configuration& par, ScalarAffordance affordance,
                           const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx =
                               std::nullopt) const -> ScalarType {
    using namespace Dune::DerivativeDirections;
    using namespace Dune;

    const auto eps     = strainFunction(par, dx);
    const auto& lambda = par.parameter();
    ScalarType energy  = 0.0;

    for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
      const auto EVoigt         = (eps.evaluate(gpIndex, on(gridElement))).eval();
      const auto internalEnergy = getInternalEnergy(EVoigt);
      energy += internalEnergy * geo_->integrationElement(gp.position()) * gp.weight();
    }

    return energy;
  }

  template <typename ScalarType>
  void calculateVectorImpl(
      const Configuration& par, VectorAffordance affordance, typename Traits::template VectorType<ScalarType> force,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto eps = strainFunction(par, dx);

    // Internal forces
    for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
      const double intElement = geo_->integrationElement(gp.position()) * gp.weight();
      const auto EVoigt       = (eps.evaluate(gpIndex, on(gridElement))).eval();
      const auto stresses     = getStress(EVoigt);
      for (size_t i = 0; i < numberOfNodes_; ++i) {
        const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
        force.template segment<myDim>(myDim * i) += bopI.transpose() * stresses * intElement;
      }
    }
  }
};

/**
 * \brief A helper function to create a non-linear elastic pre finite element.
 * \tparam MAT Type of the material.
 * \param mat Material parameters for the non-linear elastic element.
 * \return A non-linear elastic pre finite element.
 */
template <typename MAT>
auto nonLinearElastic(const MAT& mat) {
  NonLinearElasticPre<MAT> pre(mat);

  return pre;
}

} // namespace Ikarus

#else
  #error NonLinearElastic depends on dune-localfefunctions, which is not included
#endif
