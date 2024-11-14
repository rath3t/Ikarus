// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file enhancedassumedstrains.hh
 * \brief Definition of the EAS class.
 * \ingroup  mechanics
 */

#pragma once

#if HAVE_DUNE_LOCALFEFUNCTIONS
  #include <dune/localfefunctions/derivativetransformators.hh>
  #include <dune/localfefunctions/linearAlgebraHelper.hh>
  #include <dune/localfefunctions/meta.hh>

  #include <ikarus/finiteelements/fehelper.hh>
  #include <ikarus/finiteelements/ferequirements.hh>
  #include <ikarus/finiteelements/mechanics/easvariants.hh>
  #include <ikarus/finiteelements/mechanics/materials/tags.hh>
  #include <ikarus/utils/concepts.hh>

namespace Ikarus {

template <typename PreFE, typename FE, StrainTags ES>
class EnhancedAssumedStrains;

/**
 * \brief A PreFE struct for Enhanced Assumed Strains.
 * \tparam ES The strain tag that is enhanced.
 */
template <StrainTags ES>
struct EnhancedAssumedStrainsPre
{
  int numberOfParameters{};

  template <typename PreFE, typename FE>
  using Skill = EnhancedAssumedStrains<PreFE, FE, ES>;
};

/**
 * \brief Wrapper class for using Enhanced Assumed Strains (EAS) with displacement based elements.
 *
 * \ingroup mechanics
 *
 * This class extends a displacement-based element to support Enhanced Assumed Strains.
 *
 * \tparam PreFE Type of the pre finite element.
 * \tparam FE Type of the finite element.
 * \tparam ES The strain tag that is enhanced.
 */
template <typename PreFE, typename FE, StrainTags ES>
class EnhancedAssumedStrains
{
public:
  using Traits = PreFE::Traits;
  using Requirement =
      FERequirementsFactory<FESolutions::displacement, FEParameter::loadfactor, Traits::useEigenRef>::type;
  using LocalView = typename Traits::LocalView;
  using Geometry  = typename Traits::Geometry;
  using GridView  = typename Traits::GridView;
  using Pre       = EnhancedAssumedStrainsPre<ES>;

  static constexpr auto enhancedStrain = ES;

  template <typename ST>
  using VectorXOptRef = std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>;

  /**
   * \brief Constructor for Enhanced Assumed Strains elements.
   * \param pre The pre finite element
   */
  explicit EnhancedAssumedStrains(const Pre& pre) { this->setEASType(pre.numberOfParameters); }

  /**
   * \brief Checks if the element is displacement-based and the EAS is turned off.
   *
   * \return True if the element is displacement-based, false otherwise.
   */
  bool isDisplacementBased() const { return easVariant_.isDisplacmentBased(); }

  /**
   * \brief Gets the variant representing the type of Enhanced Assumed Strains (EAS).
   *
   * \return Const reference to the EAS variant.
   */
  const auto& easVariant() const { return easVariant_; }

  /**
   * \brief Gets the number of EAS parameters based on the current EAS type.
   *
   * \return Number of EAS parameters.
   */
  auto numberOfEASParameters() const { return easVariant_.numberOfEASParameters(); }

  /**
   * \brief Calculates a requested result at a specific local position using the Enhanced Assumed Strains (EAS)
   * method.
   *
   * This function calculates the results at the specified local coordinates .
   * It takes into account the displacement-based element calculations and, if applicable, incorporates the EAS method
   * for enhanced accuracy.
   *
   * \param req The result requirements.
   * \param local The local coordinates at which results are to be calculated.
   * \return calculated result
   *
   * \tparam RT The type representing the requested result.
   */
  template <template <typename, int, int> class RT>
  auto calculateAtImpl(const Requirement& req, const Dune::FieldVector<double, Traits::mydim>& local,
                       Dune::PriorityTag<2>) const {
    if (isDisplacementBased())
      return underlying().template calculateAtImpl<RT>(req, local, Dune::PriorityTag<1>());

    auto resultVector   = underlying().template calculateAtImpl<RT>(req, local, Dune::PriorityTag<1>());
    auto strainFunction = underlying().strainFunction(req);
    const auto ufunc    = underlying().displacementFunction(req);
    auto disp           = Dune::viewAsFlatEigenVector(ufunc.coefficientsRef());
    if constexpr (isSameResultType<RT, ResultTypes::linearStress>) {
      using RTWrapper = ResultWrapper<RT<typename Traits::ctype, Traits::mydim, Traits::worlddim>, ResultShape::Vector>;
      RTWrapper resultWrapper;
      auto calculateAtContribution = [&]<typename EAST>(const EAST& easFunction) {
        typename EAST::DType D;
        calculateDAndLMatrix(easFunction, req, D, L_);
        this->alpha_              = (-D.inverse() * L_ * disp).eval();
        const auto enhancedStrain = enhancedStrainFunction(easFunction, strainFunction, local);
        resultWrapper             = resultVector.asVec() + underlying().stress(enhancedStrain).eval();
      };
      easVariant_(calculateAtContribution);
      return resultWrapper;
    } else {
      DUNE_THROW(Dune::NotImplemented, "EAS method is not implemented to compute the requested result type.");
    }
  }

  /**
   * \brief Sets the EAS type for 2D elements.
   *
   * \param numberOfEASParameters The number of EAS parameters
   */
  void setEASType(int numberOfEASParameters) {
    if (numberOfEASParameters != 0)
      easApplicabilityCheck();
    easVariant_.setEASType(numberOfEASParameters);
    initializeState();
  }

  /**
   * \brief Gets the internal state variable alpha for the EAS element.
   *
   * \return Internal state variable (alpha).
   */
  const auto& alpha() const { return alpha_; }

protected:
  void bindImpl() {
    assert(underlying().localView().bound());
    easVariant_.bind(underlying().localView().element().geometry());
    initializeState();
  }

  /**
   * \brief Updates the internal state variable alpha_ at the end of an iteration before the update of the displacements
   * done by the non-linear solver. See \cite bieberLockingHourglassingNonlinear2024 for implementation details and
   * further references.
   *
   * \param par The Requirement object.
   * \param correction The correction in displacement (DeltaD) vector passed based on which the internal state variable
   * alpha is to be updated.
   */
  void updateStateImpl(const Requirement& par, typename Traits::template VectorTypeConst<> correction) const {
    using ScalarType = Traits::ctype;
    easApplicabilityCheck();
    if (isDisplacementBased())
      return;
    const auto& Rtilde      = calculateRtilde<ScalarType>(par);
    const auto localdxBlock = Ikarus::FEHelper::localSolutionBlockVector<Traits, Eigen::VectorXd, double>(
        correction, underlying().localView());
    const auto localdx  = Dune::viewAsFlatEigenVector(localdxBlock);
    decltype(auto) LMat = LMatFunc<ScalarType>();

    auto correctAlpha = [&]<typename EAST>(const EAST& easFunction) {
      constexpr int enhancedStrainSize = EAST::enhancedStrainSize;
      Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize> D;
      calculateDAndLMatrix(easFunction, par, D, LMat);
      const decltype(alpha_) updateAlpha = D.inverse() * (Rtilde + (LMat * localdx).eval());
      this->alpha_ -= updateAlpha;
    };

    easVariant_(correctAlpha);
  }

  inline void easApplicabilityCheck() const {
    const auto& numberOfNodes = underlying().numberOfNodes();
    if (not((numberOfNodes == 4 and Traits::mydim == 2) or (numberOfNodes == 8 and Traits::mydim == 3)) and
        (not isDisplacementBased()))
      DUNE_THROW(Dune::NotImplemented, "EAS is only supported for Q1 or H1 elements");
  }

  template <typename ScalarType>
  void calculateMatrixImpl(const Requirement& par, const MatrixAffordance& affordance,
                           typename Traits::template MatrixType<> K,
                           const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    if (affordance != MatrixAffordance::stiffness)
      DUNE_THROW(Dune::NotImplemented, "MatrixAffordance not implemented: " + toString(affordance));
    easApplicabilityCheck();
    if (isDisplacementBased())
      return;

    auto strainFunction  = underlying().strainFunction(par, dx);
    const auto kFunction = underlying().template stiffnessMatrixFunction<ScalarType>(par, K, dx);
    decltype(auto) LMat  = LMatFunc<ScalarType>();

    auto calculateMatrixContribution = [&]<typename EAST>(const EAST& easFunction) {
      typename EAST::DType D;
      calculateDAndLMatrix(easFunction, par, D, LMat, dx);

      for (const auto& [gpIndex, gp] : strainFunction.viewOverIntegrationPoints()) {
        const auto enhancedStrain = enhancedStrainFunction(easFunction, strainFunction, gp.position());
        kFunction(enhancedStrain, gpIndex, gp);
      }

      K.template triangularView<Eigen::Upper>() -= LMat.transpose() * D.inverse() * LMat;
      K.template triangularView<Eigen::StrictlyLower>() = K.transpose();
    };
    easVariant_(calculateMatrixContribution);
  }

  template <typename ScalarType>
  inline ScalarType calculateScalarImpl(const Requirement& par, ScalarAffordance affordance,
                                        const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    easApplicabilityCheck();
    if (isDisplacementBased())
      return 0.0;
    DUNE_THROW(Dune::NotImplemented,
               "EAS element do not support any scalar calculations, i.e. they are not derivable from a potential");
  }

  template <typename ScalarType>
  void calculateVectorImpl(const Requirement& par, VectorAffordance affordance,
                           typename Traits::template VectorType<ScalarType> force,
                           const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    if (affordance != VectorAffordance::forces)
      DUNE_THROW(Dune::NotImplemented, "VectorAffordance not implemented: " + toString(affordance));
    easApplicabilityCheck();
    if (isDisplacementBased())
      return;
    auto strainFunction     = underlying().strainFunction(par, dx);
    const auto fIntFunction = underlying().template internalForcesFunction<ScalarType>(par, force, dx);

    auto calculateForceContribution = [&]<typename EAST>(const EAST& easFunction) {
      typename EAST::DType D;
      calculateDAndLMatrix(easFunction, par, D, L_);

      const auto& Rtilde = calculateRtilde(par, dx);
      for (const auto& [gpIndex, gp] : strainFunction.viewOverIntegrationPoints()) {
        const auto enhancedStrain = enhancedStrainFunction(easFunction, strainFunction, gp.position());
        fIntFunction(enhancedStrain, gpIndex, gp);
      }
      force -= L_.transpose() * D.inverse() * Rtilde;
    };
    easVariant_(calculateForceContribution);
  }

private:
  EAS::Impl::EASVariant<Geometry> easVariant_;
  mutable Eigen::MatrixXd L_;
  mutable Eigen::VectorXd alpha_;

  //> CRTP
  const auto& underlying() const { return static_cast<const FE&>(*this); }
  auto& underlying() { return static_cast<FE&>(*this); }

  /**
   * \brief Initializes the internal state variable alpha_ based on the number of EAS parameters.
   */
  void initializeState() {
    if (isDisplacementBased())
      return;
    alpha_.resize(numberOfEASParameters());
    alpha_.setZero();
  }

  template <typename ScalarType>
  decltype(auto) LMatFunc() const {
    if constexpr (Concepts::AutodiffScalar<ScalarType>)
      return Eigen::MatrixX<ScalarType>{};
    else
      return [this]() -> Eigen::MatrixXd& { return L_; }();
  }

  template <typename EAST, typename SF, typename POS>
  auto enhancedStrainFunction(const EAST& easFunction, const SF& strainFunction, const POS& gpPos) const {
    if constexpr (Concepts::Formulations::TotalLagrangian<FE::strainType, FE::stressType>) {
      if constexpr (enhancedStrain == StrainTags::linear or enhancedStrain == StrainTags::greenLagrangian) {
        using namespace Dune;
        using namespace Dune::DerivativeDirections;
        const auto M              = easFunction.calcM(gpPos);
        const auto EVoigt         = strainFunction.evaluate(gpPos, on(gridElement));
        const auto easStrain      = (M * alpha_).eval();
        const auto enhancedStrain = (EVoigt + easStrain).eval();
        return enhancedStrain;
      } else {
        static_assert(Dune::AlwaysFalse<FE>::value,
                      "EAS method is not implemented for the provided strain tag that is to be enhanced.");
      }
    } else {
      static_assert(Dune::AlwaysFalse<FE>::value,
                    "EAS method is not implemented for the provided stress and strain measure.");
    }
  }

  template <typename ScalarType, int enhancedStrainSize>
  void calculateDAndLMatrix(const auto& easFunction, const auto& par,
                            Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize>& DMat,
                            Eigen::MatrixX<ScalarType>& LMat,
                            const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    using namespace Dune;
    using namespace Dune::DerivativeDirections;

    auto strainFunction      = underlying().strainFunction(par, dx);
    const auto geo           = underlying().localView().element().geometry();
    const auto numberOfNodes = underlying().numberOfNodes();
    DMat.setZero();
    LMat.setZero(enhancedStrainSize, underlying().localView().size());
    for (const auto& [gpIndex, gp] : strainFunction.viewOverIntegrationPoints()) {
      const auto enhancedStrain = enhancedStrainFunction(easFunction, strainFunction, gp.position());
      const auto M              = easFunction.calcM(gp.position());
      const auto CEval          = underlying().materialTangent(enhancedStrain);
      const double intElement   = geo.integrationElement(gp.position()) * gp.weight();
      DMat += M.transpose() * CEval * M * intElement;
      for (size_t i = 0U; i < numberOfNodes; ++i) {
        const size_t I = Traits::worlddim * i;
        const auto Bi  = strainFunction.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
        LMat.template block<enhancedStrainSize, Traits::worlddim>(0, I) += M.transpose() * CEval * Bi * intElement;
      }
    }
  }

  template <typename ScalarType>
  Eigen::VectorX<ScalarType> calculateRtilde(const Requirement& par,
                                             const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    const auto geo      = underlying().localView().element().geometry();
    auto strainFunction = underlying().strainFunction(par, dx);
    Eigen::VectorX<ScalarType> Rtilde;
    Rtilde.setZero(numberOfEASParameters());

    auto calculateRtildeContribution = [&]<typename EAST>(const EAST& easFunction) {
      for (const auto& [gpIndex, gp] : strainFunction.viewOverIntegrationPoints()) {
        const auto enhancedStrain = enhancedStrainFunction(easFunction, strainFunction, gp.position());
        const auto M              = easFunction.calcM(gp.position());
        const double intElement   = geo.integrationElement(gp.position()) * gp.weight();
        auto stresses             = underlying().stress(enhancedStrain);
        Rtilde += (M.transpose() * stresses).eval() * intElement;
      }
    };

    easVariant_(calculateRtildeContribution);
    return Rtilde;
  }
};

/**
 * \brief A helper function to create an enhanced assumed strain pre finite element.
 * \tparam ES The strain tag that is enhanced.
 * \param numberOfEASParameters Number of EAS parameters
 * \return An enhanced assumed strain pre finite element.
 */
template <StrainTags ES = StrainTags::linear>
auto eas(int numberOfEASParameters = 0) {
  EnhancedAssumedStrainsPre<ES> pre(numberOfEASParameters);

  return pre;
}

} // namespace Ikarus

#else
  #error EnhancedAssumedStrains depends on dune-localfefunctions, which is not included
#endif
