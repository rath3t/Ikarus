// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file kirchhoffloveshell.hh
 * \brief Definition of the KirchhoffLoveShell class for Kirchhoff-Love shell elements in Ikarus.
 */

#pragma once

#include <dune/geometry/quadraturerules.hh>
#include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
#include <dune/localfefunctions/impl/standardLocalFunction.hh>

#include <ikarus/finiteelements/fehelper.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/feresulttypes.hh>
#include <ikarus/finiteelements/mechanics/loads.hh>
#include <ikarus/finiteelements/mechanics/membranestrains.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/linearalgebrahelper.hh>

namespace Ikarus {
template <typename PreFE, typename FE>
class KirchhoffLoveShell;

/**
 * \brief A PreFE struct for Kirchhoff-Love shell elements.
 */
struct KirchhoffLoveShellPre
{
  YoungsModulusAndPoissonsRatio material;
  double thickness;

  template <typename PreFE, typename FE>
  using Skill = KirchhoffLoveShell<PreFE, FE>;
};

/**
 * \brief Kirchhoff-Love shell finite element class.
 *
 * This class represents Kirchhoff-Love shell finite elements.
 *
 * \ingroup mechanics
 *
 * \tparam PreFE Type of the pre finite element.
 * \tparam FE Type of the finite element.
 */
template <typename PreFE, typename FE>
class KirchhoffLoveShell : public ResultTypeBase<>
{
public:
  using Traits       = PreFE::Traits;
  using BasisHandler = typename Traits::BasisHandler;
  using FlatBasis    = typename Traits::FlatBasis;
  using Requirement =
      FERequirementsFactory<FESolutions::displacement, FEParameter::loadfactor, Traits::useEigenRef>::type;
  using LocalView      = typename Traits::LocalView;
  using Geometry       = typename Traits::Geometry;
  using GridView       = typename Traits::GridView;
  using Element        = typename Traits::Element;
  using LocalBasisType = decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis());
  using Pre            = KirchhoffLoveShellPre;

  static constexpr int myDim              = Traits::mydim;
  static constexpr int worldDim           = Traits::worlddim;
  static constexpr int membraneStrainSize = 3;
  static constexpr int bendingStrainSize  = 3;

  /**
   * \brief A structure representing kinematic variables.
   *
   * This structure holds various kinematic variables used in a mechanical analysis. It includes material tangent,
   * membrane strain, bending strain, Jacobian matrices of deformed and reference geometries, Hessian matrices of
   * deformed and reference geometries, the normal vector, and the normalized normal vector.
   *
   * \tparam ST The scalar type for the matrix and vector elements.
   */
  template <typename ST = double>
  struct KinematicVariables
  {
    Eigen::Matrix<double, 3, 3> C; ///< material tangent
    Eigen::Vector3<ST> epsV;       ///< membrane strain in Voigt notation
    Eigen::Vector3<ST> kappaV;     ///< bending strain in Voigt notation
    Eigen::Matrix<ST, 2, 3> j;     ///< Jacobian of the deformed geometry
    Eigen::Matrix<double, 2, 3> J; ///< Jacobian of the reference geometry
    Eigen::Matrix3<ST> h;          ///< Hessian of the deformed geometry
    Eigen::Matrix3<double> H;      ///< Hessian of the reference geometry
    Eigen::Vector3<ST> a3N;        ///< Normal vector of the deformed geometry
    Eigen::Vector3<ST> a3;         ///< normalized normal vector of the deformed geometry
  };

  /**
   * \brief Constructor for the KirchhoffLoveShell class.
   *
   * Initializes the KirchhoffLoveShell instance with the given parameters.
   *
   * \param pre The pre fe
   */
  KirchhoffLoveShell(const Pre& pre)
      : mat_{pre.material},
        thickness_{pre.thickness} {}

protected:
  /**
   * \brief A helper function to bind the local view to the element.
   */
  void bindImpl() {
    const auto& localView = underlying().localView();
    assert(localView.bound());
    const auto& element = localView.element();
    auto& firstChild    = localView.tree().child(0);

    const auto& fe = firstChild.finiteElement();
    // geo_.emplace(element.geometry());
    numberOfNodes_ = fe.size();
    order_         = 2 * (fe.localBasis().order());
    localBasis_    = Dune::CachedLocalBasis(fe.localBasis());
    if constexpr (requires { element.impl().getQuadratureRule(order_); })
      if (element.impl().isTrimmed())
        localBasis_.bind(element.impl().getQuadratureRule(order_), Dune::bindDerivatives(0, 1, 2));
      else
        localBasis_.bind(Dune::QuadratureRules<double, myDim>::rule(element.type(), order_),
                         Dune::bindDerivatives(0, 1, 2));
    else
      localBasis_.bind(Dune::QuadratureRules<double, myDim>::rule(element.type(), order_),
                       Dune::bindDerivatives(0, 1, 2));
  }

public:
  /**
   * \brief Get the displacement function and nodal displacements.
   *
   * Retrieves the displacement function and nodal displacements based on the given FERequirements.
   *
   * \tparam ST The scalar type used for calculations.
   * \param par The FERequirements.
   * \param dx Optional additional displacement vector.
   * \return The displacement function.
   */
  template <typename ST = double>
  auto displacementFunction(
      const Requirement& par,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>& dx = std::nullopt) const {
    const auto& d = par.globalSolution();
    auto disp     = Ikarus::FEHelper::localSolutionBlockVector<Traits>(d, underlying().localView(), dx);
    Dune::StandardLocalFunction uFunction(
        localBasis_, disp, std::make_shared<const Geometry>(underlying().localView().element().geometry()));
    return uFunction;
  }

  Geometry geometry() const { return underlying().localView().element().geometry(); }
  [[nodiscard]] size_t numberOfNodes() const { return numberOfNodes_; }
  [[nodiscard]] int order() const { return order_; }

  /**
   * \brief Calculates a requested result at a specific local position.
   *
   * \param req The FERequirementType object holding the global solution.
   * \param local Local position vector.
   * \return calculated result
   *
   * \tparam RT The type representing the requested result.
   */
  template <template <typename, int, int> class RT>
  requires(supportsResultType<RT>())
  auto calculateAtImpl([[maybe_unused]] const Requirement& req,
                       [[maybe_unused]] const Dune::FieldVector<double, Traits::mydim>& local)
      -> ResultWrapper<RT<double, myDim, worldDim>, ResultShape::Vector> {
    DUNE_THROW(Dune::NotImplemented, "No results are implemented");
  }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const FE&>(*this); }
  auto& underlying() { return static_cast<FE&>(*this); }
  // std::optional<Geometry> geo_;
  Dune::CachedLocalBasis<std::remove_cvref_t<LocalBasisType>> localBasis_;
  // DefaultMembraneStrain membraneStrain_;
  YoungsModulusAndPoissonsRatio mat_;
  double thickness_;

  size_t numberOfNodes_{0};
  int order_{};

protected:
  /**
   * \brief Compute material properties and strains at a given integration point.
   *
   * \param gpPos The position of the integration point.
   * \param gpIndex The index of the integration point.
   * \param geo The geometry object providing position and derivatives.
   * \param uFunction The function representing the displacement field.
   *
   * \return A tuple containing the material tangent, membrane strain, bending,
   *         Jacobian matrix of the reference position,  Jacobian matrix of the current position, Hessian matrix of
   * the current position, Hessian matrix of
   * the reference position, normal vector, and normalized normal vector at the given
   * integration point.
   */
  auto computeMaterialAndStrains(const Dune::FieldVector<double, 2>& gpPos, int gpIndex, const Geometry& geo,
                                 const auto& uFunction) const {
    using ST = typename std::remove_cvref_t<decltype(uFunction)>::ctype;

    KinematicVariables<ST> kin;
    using namespace Dune;
    using namespace Dune::DerivativeDirections;
    const auto [X, Jd, Hd]              = geo.impl().zeroFirstAndSecondDerivativeOfPosition(gpPos);
    kin.J                               = toEigen(Jd);
    kin.H                               = toEigen(Hd);
    const Eigen::Matrix<double, 2, 2> A = kin.J * kin.J.transpose();
    Eigen::Matrix<double, 3, 3> G       = Eigen::Matrix<double, 3, 3>::Zero();

    G.block<2, 2>(0, 0)                    = A;
    G(2, 2)                                = 1;
    const Eigen::Matrix<double, 3, 3> GInv = G.inverse();
    kin.C                                  = materialTangent(GInv);

    kin.epsV = DefaultMembraneStrain::value(gpPos, geo, uFunction);

    const auto& Ndd                     = localBasis_.evaluateSecondDerivatives(gpIndex);
    const auto uasMatrix                = Dune::viewAsEigenMatrixAsDynFixed(uFunction.coefficientsRef());
    const auto hessianu                 = Ndd.transpose().template cast<ST>() * uasMatrix;
    kin.h                               = kin.H + hessianu;
    const Eigen::Matrix<ST, 3, 2> gradu = toEigen(uFunction.evaluateDerivative(
        gpIndex, Dune::wrt(spatialAll), Dune::on(Dune::DerivativeDirections::referenceElement)));
    kin.j                               = kin.J + gradu.transpose();
    kin.a3N                             = (kin.j.row(0).cross(kin.j.row(1)));
    kin.a3                              = kin.a3N.normalized();
    Eigen::Vector<ST, 3> bV             = kin.h * kin.a3;
    bV(2) *= 2; // Voigt notation requires the two here
    const auto BV = toVoigt(toEigen(geo.impl().secondFundamentalForm(gpPos)));
    kin.kappaV    = BV - bV;
    return kin;
  }

  void updateStateImpl(const Requirement& /* par */,
                       typename Traits::template VectorTypeConst<> /* correction */) const {}

  template <typename ST>
  void calculateMatrixImpl(
      const Requirement& par, const MatrixAffordance& affordance, typename Traits::template MatrixType<ST> K,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>& dx = std::nullopt) const {
    if (affordance != MatrixAffordance::stiffness)
      DUNE_THROW(Dune::NotImplemented, "MatrixAffordance not implemented: " + toString(affordance));
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto uFunction = displacementFunction(par, dx);
    const auto& lambda   = par.parameter();
    const auto geo       = underlying().localView().element().geometry();

    for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
      const auto intElement = geo.integrationElement(gp.position()) * gp.weight();
      const auto [C, epsV, kappaV, jE, J, h, H, a3N, a3] =
          computeMaterialAndStrains(gp.position(), gpIndex, geo, uFunction);
      const Eigen::Vector<ST, membraneStrainSize> membraneForces = thickness_ * C * epsV;
      const Eigen::Vector<ST, bendingStrainSize> moments         = Dune::power(thickness_, 3) / 12.0 * C * kappaV;

      const auto& Nd  = localBasis_.evaluateJacobian(gpIndex);
      const auto& Ndd = localBasis_.evaluateSecondDerivatives(gpIndex);
      for (size_t i = 0; i < numberOfNodes_; ++i) {
        Eigen::Matrix<ST, membraneStrainSize, worldDim> bopIMembrane =
            DefaultMembraneStrain::derivative(gp.position(), jE, Nd, geo, uFunction, localBasis_, i);

        Eigen::Matrix<ST, bendingStrainSize, worldDim> bopIBending = bopBending(jE, h, Nd, Ndd, i, a3N, a3);
        for (size_t j = i; j < numberOfNodes_; ++j) {
          auto KBlock = K.template block<worldDim, worldDim>(worldDim * i, worldDim * j);
          Eigen::Matrix<ST, membraneStrainSize, worldDim> bopJMembrane =
              DefaultMembraneStrain::derivative(gp.position(), jE, Nd, geo, uFunction, localBasis_, j);
          Eigen::Matrix<ST, bendingStrainSize, worldDim> bopJBending = bopBending(jE, h, Nd, Ndd, j, a3N, a3);
          KBlock += thickness_ * bopIMembrane.transpose() * C * bopJMembrane * intElement;
          KBlock += Dune::power(thickness_, 3) / 12.0 * bopIBending.transpose() * C * bopJBending * intElement;

          Eigen::Matrix<ST, worldDim, worldDim> kgMembraneIJ = DefaultMembraneStrain::secondDerivative(
              gp.position(), Nd, geo, uFunction, localBasis_, membraneForces, i, j);
          Eigen::Matrix<ST, worldDim, worldDim> kgBendingIJ = kgBending(jE, h, Nd, Ndd, a3N, a3, moments, i, j);
          KBlock += kgMembraneIJ * intElement;
          KBlock += kgBendingIJ * intElement;
        }
      }
    }
    K.template triangularView<Eigen::StrictlyLower>() = K.transpose();
  }

  template <typename ST>
  void calculateVectorImpl(
      const Requirement& par, const VectorAffordance& affordance, typename Traits::template VectorType<ST> force,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>& dx = std::nullopt) const {
    if (affordance != VectorAffordance::forces)
      DUNE_THROW(Dune::NotImplemented, "VectorAffordance not implemented: " + toString(affordance));
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto uFunction = displacementFunction(par, dx);
    const auto& lambda   = par.parameter();
    const auto geo       = underlying().localView().element().geometry();

    // Internal forces
    for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
      const auto [C, epsV, kappaV, jE, J, h, H, a3N, a3] =
          computeMaterialAndStrains(gp.position(), gpIndex, geo, uFunction);
      const Eigen::Vector<ST, 3> membraneForces = thickness_ * C * epsV;
      const Eigen::Vector<ST, 3> moments        = Dune::power(thickness_, 3) / 12.0 * C * kappaV;

      const auto& Nd  = localBasis_.evaluateJacobian(gpIndex);
      const auto& Ndd = localBasis_.evaluateSecondDerivatives(gpIndex);
      for (size_t i = 0; i < numberOfNodes_; ++i) {
        Eigen::Matrix<ST, 3, 3> bopIMembrane =
            DefaultMembraneStrain::derivative(gp.position(), jE, Nd, geo, uFunction, localBasis_, i);
        Eigen::Matrix<ST, 3, 3> bopIBending = bopBending(jE, h, Nd, Ndd, i, a3N, a3);
        force.template segment<3>(3 * i) +=
            bopIMembrane.transpose() * membraneForces * geo.integrationElement(gp.position()) * gp.weight();
        force.template segment<3>(3 * i) +=
            bopIBending.transpose() * moments * geo.integrationElement(gp.position()) * gp.weight();
      }
    }
  }

  template <typename ST>
  auto calculateScalarImpl(
      const Requirement& par, const ScalarAffordance& affordance,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>& dx = std::nullopt) const -> ST {
    if (affordance != ScalarAffordance::mechanicalPotentialEnergy)
      DUNE_THROW(Dune::NotImplemented, "ScalarAffordance not implemented: " + toString(affordance));
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto uFunction = displacementFunction(par, dx);
    const auto& lambda   = par.parameter();
    ST energy            = 0.0;

    const auto geo = underlying().localView().element().geometry();

    for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
      const auto [C, epsV, kappaV, j, J, h, H, a3N, a3] =
          computeMaterialAndStrains(gp.position(), gpIndex, geo, uFunction);

      const ST membraneEnergy = 0.5 * thickness_ * epsV.dot(C * epsV);
      const ST bendingEnergy  = 0.5 * Dune::power(thickness_, 3) / 12.0 * kappaV.dot(C * kappaV);
      energy += (membraneEnergy + bendingEnergy) * geo.integrationElement(gp.position()) * gp.weight();
    }

    return energy;
  }

  template <typename ST>
  Eigen::Matrix<ST, 3, 3> kgBending(const Eigen::Matrix<ST, 2, 3>& jcur, const Eigen::Matrix3<ST>& h, const auto& dN,
                                    const auto& ddN, const Eigen::Vector3<ST>& a3N, const Eigen::Vector3<ST>& a3,
                                    const Eigen::Vector3<ST>& S, int I, int J) const {
    Eigen::Matrix<ST, 3, 3> kg;
    kg.setZero();

    const auto& dN1i = dN(I, 0);
    const auto& dN1j = dN(J, 0);
    const auto& dN2i = dN(I, 1);
    const auto& dN2j = dN(J, 1);

    const Eigen::Matrix<ST, 3, 3> P =
        1.0 / a3N.norm() * (Eigen::Matrix<double, 3, 3>::Identity() - a3 * a3.transpose());

    const auto a1dxI = Eigen::Matrix<double, 3, 3>::Identity() *
                       dN1i; // the auto here allows the exploitation of the identity matrices,
    // due to Eigen's expression templates
    const auto a2dxI                    = Eigen::Matrix<double, 3, 3>::Identity() * dN2i;
    const auto a1dxJ                    = Eigen::Matrix<double, 3, 3>::Identity() * dN1j;
    const auto a2dxJ                    = Eigen::Matrix<double, 3, 3>::Identity() * dN2j;
    const auto a1                       = jcur.row(0);
    const auto a2                       = jcur.row(1);
    const Eigen::Matrix<ST, 3, 3> a3NdI = a1dxI.colwise().cross(a2) - a2dxI.colwise().cross(a1);
    const Eigen::Matrix<ST, 3, 3> a3NdJ = a1dxJ.colwise().cross(a2) - a2dxJ.colwise().cross(a1);
    Eigen::Matrix<ST, 3, 3> a3dI        = P * a3NdI;
    Eigen::Matrix<ST, 3, 3> a3dJ        = P * a3NdJ;
    for (int i = 0; i < 3; ++i) {
      const auto a_albe       = h.row(i).transpose();
      const auto& ddNI        = ddN(I, i);
      const auto& ddNJ        = ddN(J, i);
      Eigen::Vector3<ST> vecd = P * a_albe;

      Eigen::Matrix<ST, 3, 3> a3Ndd =
          1.0 / a3N.squaredNorm() *
          ((3 * a3 * a3.transpose() - Eigen::Matrix<double, 3, 3>::Identity()) * (a3.dot(a_albe)) -
           a_albe * a3.transpose() - a3 * a_albe.transpose());

      Eigen::Matrix<ST, 3, 3> secondDerivativeDirectorIJ = skew(((dN2i * dN1j - dN1i * dN2j) * vecd).eval());
      kg -= (a3NdI.transpose() * a3Ndd * a3NdJ + secondDerivativeDirectorIJ + (ddNI * a3dJ + ddNJ * a3dI.transpose())) *
            S[i] * (i == 2 ? 2 : 1);
    }

    return kg;
  }

  template <typename ST>
  Eigen::Matrix<ST, 3, 3> bopBending(const Eigen::Matrix<ST, 2, 3>& jcur, const Eigen::Matrix3<ST>& h, const auto& dN,
                                     const auto& ddN, const int node, const Eigen::Vector3<ST>& a3N,
                                     const Eigen::Vector3<ST>& a3) const {
    const Eigen::Matrix<ST, 3, 3> a1dxI =
        Eigen::Matrix<double, 3, 3>::Identity() * dN(node, 0); // this should be double
    // but the cross-product below complains otherwise
    const Eigen::Matrix<ST, 3, 3> a2dxI = Eigen::Matrix<double, 3, 3>::Identity() * dN(node, 1);
    const auto a1                       = jcur.row(0);
    const auto a2                       = jcur.row(1);
    const Eigen::Matrix<ST, 3, 3> a3NdI =
        a1dxI.colwise().cross(a2) - a2dxI.colwise().cross(a1); // minus needed since order has
    // to be swapped to get column-wise cross product working
    const Eigen::Matrix<ST, 3, 3> a3d1 =
        1.0 / a3N.norm() * (Eigen::Matrix<double, 3, 3>::Identity() - a3 * a3.transpose()) * a3NdI;

    Eigen::Matrix<ST, 3, 3> bop = -(h * a3d1 + (a3 * ddN.row(node)).transpose());
    bop.row(2) *= 2;

    return bop;
  }

  /**
   * \brief Gets the material tangent matrix for the linear elastic material.
   * \param Aconv A tensor to get the material tangent in the curvilinear coordinate system.
   *
   * \return The material tangent matrix.
   */
  Eigen::Matrix<double, 3, 3> materialTangent(const Eigen::Matrix<double, 3, 3>& Aconv) const {
    const double lambda   = mat_.emodul * mat_.nu / ((1.0 + mat_.nu) * (1.0 - 2.0 * mat_.nu));
    const double mu       = mat_.emodul / (2.0 * (1.0 + mat_.nu));
    const double lambdbar = 2.0 * lambda * mu / (lambda + 2.0 * mu);
    Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3, 3>> moduli;
    const auto AconvT = tensorView(Aconv, std::array<Eigen::Index, 2>({3, 3}));
    moduli = lambdbar * dyadic(AconvT, AconvT).eval() + 2.0 * mu * symmetricFourthOrder<double>(Aconv, Aconv);

    auto C                          = toVoigt(moduli);
    Eigen::Matrix<double, 3, 3> C33 = C({0, 1, 5}, {0, 1, 5});
    return C33;
  }
};

/**
 * \brief A struct containing information about the Youngs Modulus,
 * Poisson's ratio and the thickness for the Kirchhoff-Love shell element.
 */
struct KlArgs
{
  double youngs_modulus;
  double nu;
  double thickness;
};

/**
 * \brief A helper function to create a Kirchhoff-Love shell pre finite element.
 * \param args Arguments for the Kirchhoff-Love shell element.
 * \return A Kirchhoff-Love shell pre finite element.
 */
auto kirchhoffLoveShell(const KlArgs& args) {
  KirchhoffLoveShellPre pre({.emodul = args.youngs_modulus, .nu = args.nu}, args.thickness);

  return pre;
}
} // namespace Ikarus
