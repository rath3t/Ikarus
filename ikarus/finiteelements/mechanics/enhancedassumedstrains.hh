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
  #include <ikarus/finiteelements/mixin.hh>
  #include <ikarus/utils/concepts.hh>

namespace Ikarus {

template <typename PreFE, typename FE>
class EnhancedAssumedStrains;

/**
 * \brief A PreFE struct for Enhanced Assumed Strains.
 */
struct EnhancedAssumedStrainsPre
{
  int numberOfParameters{};

  template <typename PreFE, typename FE>
  using Skill = EnhancedAssumedStrains<PreFE, FE>;
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
 */
template <typename PreFE, typename FE>
class EnhancedAssumedStrains
{
public:
  using Traits = PreFE::Traits;
  using Requirement =
      FERequirementsFactory<FESolutions::displacement, FEParameter::loadfactor, Traits::useEigenRef>::type;
  using LocalView = typename Traits::LocalView;
  using Geometry  = typename Traits::Geometry;
  using GridView  = typename Traits::GridView;
  using Pre       = EnhancedAssumedStrainsPre;

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

    using namespace Dune::Indices;
    using namespace Dune::DerivativeDirections;
    auto resultVector = underlying().template calculateAtImpl<RT>(req, local, Dune::PriorityTag<1>());
    if constexpr (isSameResultType<RT, ResultTypes::linearStress>) {
      using RTWrapper = ResultWrapper<RT<typename Traits::ctype, Traits::mydim, Traits::worlddim>, ResultShape::Vector>;
      RTWrapper resultWrapper;

      auto calculateAtContribution = [&]<typename EAST>(const EAST& easFunction) {
        typename EAST::DType D;
        calculateDAndLMatrix(easFunction, req, D, L_);
        const auto ufunc = underlying().displacementFunction(req);
        const auto disp  = Dune::viewAsFlatEigenVector(ufunc.coefficientsRef());
        const auto C     = underlying().materialTangentFunction(req);
        const auto alpha = (-D.inverse() * L_ * disp).eval();
        const auto M     = easFunction.calcM(local);
        const auto CEval = C(local);
        auto easStress   = (CEval * M * alpha).eval();
        resultWrapper    = resultVector.asVec() + easStress;
      };
      easVariant_(calculateAtContribution);
      return resultWrapper;
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
  }

protected:
  void bindImpl() {
    assert(underlying().localView().bound());
    easVariant_.bind(underlying().localView().element().geometry());
  }

public:
protected:
  inline void easApplicabilityCheck() const {
    const auto& numberOfNodes = underlying().numberOfNodes();
    assert(not(not((numberOfNodes == 4 and Traits::mydim == 2) or (numberOfNodes == 8 and Traits::mydim == 3)) and
               (not isDisplacementBased())) &&
           "EAS only supported for Q1 or H1 elements");
  }

  template <typename ScalarType>
  void calculateMatrixImpl(
      const Requirement& par, const MatrixAffordance& affordance, typename Traits::template MatrixType<> K,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    easApplicabilityCheck();
    if (isDisplacementBased() or affordance != MatrixAffordance::stiffness)
      return;

    decltype(auto) LMat = [this]() -> decltype(auto) {
      if constexpr (Concepts::AutodiffScalar<ScalarType>)
        return Eigen::MatrixX<ScalarType>{};
      else
        return [this]() -> Eigen::MatrixXd& { return L_; }();
    }();

    auto calculateMatrixContribution = [&]<typename EAST>(const EAST& easFunction) {
      typename EAST::DType D;
      calculateDAndLMatrix(easFunction, par, D, LMat, dx);

      K.template triangularView<Eigen::Upper>() -= LMat.transpose() * D.inverse() * LMat;
      K.template triangularView<Eigen::StrictlyLower>() = K.transpose();
    };
    easVariant_(calculateMatrixContribution);
  }

  template <typename ScalarType>
  inline ScalarType calculateScalarImpl(
      const Requirement& par, ScalarAffordance affordance,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    easApplicabilityCheck();
    if (isDisplacementBased())
      return 0.0;
    DUNE_THROW(Dune::NotImplemented,
               "EAS element do not support any scalar calculations, i.e. they are not derivable from a potential");
  }

  template <typename ScalarType>
  void calculateVectorImpl(
      const Requirement& par, VectorAffordance affordance, typename Traits::template VectorType<ScalarType> force,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    easApplicabilityCheck();
    using namespace Dune;
    using namespace Dune::DerivativeDirections;
    const auto uFunction      = underlying().displacementFunction(par, dx);
    auto strainFunction       = underlying().strainFunction(par, dx);
    const auto& numberOfNodes = underlying().numberOfNodes();

    auto calculateForceContribution = [&]<typename EAST>(const EAST& easFunction) {
      typename EAST::DType D;
      calculateDAndLMatrix(easFunction, par, D, L_);

      decltype(auto) LMat = [this]() -> decltype(auto) {
        if constexpr (Concepts::AutodiffScalar<ScalarType>)
          return Eigen::MatrixX<ScalarType>{};
        else
          return [this]() -> Eigen::MatrixXd& { return L_; }();
      }();
      const auto disp  = Dune::viewAsFlatEigenVector(uFunction.coefficientsRef());
      const auto alpha = (-D.inverse() * L_ * disp).eval();
      const auto geo   = underlying().localView().element().geometry();
      auto C           = underlying().materialTangentFunction(par);

      for (const auto& [gpIndex, gp] : strainFunction.viewOverIntegrationPoints()) {
        const auto M            = easFunction.calcM(gp.position());
        const double intElement = geo.integrationElement(gp.position()) * gp.weight();
        const auto CEval        = C(gpIndex);
        auto stresses           = (CEval * M * alpha).eval();
        for (size_t i = 0; i < numberOfNodes; ++i) {
          const auto bopI = strainFunction.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
          force.template segment<Traits::worlddim>(Traits::worlddim * i) += bopI.transpose() * stresses * intElement;
        }
      }
    };
    easVariant_(calculateForceContribution);
  }

private:
  EAS::Impl::EASVariant<Geometry> easVariant_;
  mutable Eigen::MatrixXd L_;

  //> CRTP
  const auto& underlying() const { return static_cast<const FE&>(*this); }
  auto& underlying() { return static_cast<FE&>(*this); }

  template <typename ScalarType, int enhancedStrainSize>
  void calculateDAndLMatrix(
      const auto& easFunction, const auto& par, Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize>& DMat,
      Eigen::MatrixX<ScalarType>& LMat,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    using namespace Dune;
    using namespace Dune::DerivativeDirections;

    auto strainFunction      = underlying().strainFunction(par, dx);
    const auto C             = underlying().materialTangentFunction(par);
    const auto geo           = underlying().localView().element().geometry();
    const auto numberOfNodes = underlying().numberOfNodes();
    DMat.setZero();
    LMat.setZero(enhancedStrainSize, underlying().localView().size());
    for (const auto& [gpIndex, gp] : strainFunction.viewOverIntegrationPoints()) {
      const auto M            = easFunction.calcM(gp.position());
      const auto CEval        = C(gpIndex);
      const double detJTimesW = geo.integrationElement(gp.position()) * gp.weight();
      DMat += M.transpose() * CEval * M * detJTimesW;
      for (size_t i = 0U; i < numberOfNodes; ++i) {
        const size_t I = Traits::worlddim * i;
        const auto Bi  = strainFunction.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
        LMat.template block<enhancedStrainSize, Traits::worlddim>(0, I) += M.transpose() * CEval * Bi * detJTimesW;
      }
    }
  }
};

/**
 * \brief A helper function to create an enhanced assumed strain pre finite element.
 * \param numberOfEASParameters Number of EAS parameters
 * \return An enhanced assumed strain pre finite element.
 */
auto eas(int numberOfEASParameters = 0) {
  EnhancedAssumedStrainsPre pre(numberOfEASParameters);

  return pre;
}

} // namespace Ikarus

#else
  #error EnhancedAssumedStrains depends on dune-localfefunctions, which is not included
#endif
