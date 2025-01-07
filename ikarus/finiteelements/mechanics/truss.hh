// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file truss.hh
 * \brief Definition of the Truss class for finite element mechanics computations.
 * \ingroup  mechanics
 */

#pragma once

#include <optional>

#include <dune/localfefunctions/eigenDuneTransformations.hh>

#include <ikarus/finiteelements/fehelper.hh>
#include <ikarus/finiteelements/ferequirements.hh>

namespace Ikarus {

template <typename PreFE, typename FE>
class Truss;

/**
 * \brief A PreFE struct for truss elements.
 */
struct TrussPre
{
  double E;       // Young's modulus
  double A;       // Cross-section area
  double density; // Density

  template <typename PreFE, typename FE>
  using Skill = Truss<PreFE, FE>;
};

/**
 * \brief Truss class represents a truss finite element.
 *
 * \ingroup mechanics
 *
 * \tparam PreFE The type of the  pre finite element.
 * \tparam FE The type of the finite element.
 */
template <typename PreFE, typename FE>
class Truss
    : public ResultTypeBase<ResultTypes::cauchyAxialForce, ResultTypes::PK2AxialForce, ResultTypes::linearAxialForce>
{
public:
  using Traits       = PreFE::Traits;
  using BasisHandler = typename Traits::BasisHandler;
  using FlatBasis    = typename Traits::FlatBasis;
  using Requirement =
      FERequirementsFactory<FESolutions::displacement, FEParameter::loadfactor, Traits::useEigenRef>::type;
  using LocalView = typename Traits::LocalView;
  using Geometry  = typename Traits::Geometry;
  using GridView  = typename Traits::GridView;
  using Element   = typename Traits::Element;
  using Pre       = TrussPre;

  static constexpr int myDim    = Traits::mydim;
  static constexpr int worldDim = Traits::worlddim;
  static_assert(myDim == 1, "Truss elements should have myDim == 1");

  /**
   * \brief A structure representing kinematic variables.
   *
   * It includes Green-Lagrange strain, its first and second derivatives w.r.t displacements and
   * the length of the undeformed and deformed geometry.
   *
   * \tparam ST The scalar type for the matrix and vector elements.
   */
  template <typename ST>
  struct KinematicVariables
  {
    double L;                  ///< Length of the reference geometry
    Eigen::VectorXd A1;        ///< Length of the reference geometry as vector
    ST l;                      ///< Length of the deformed geometry
    ST Elin;                   ///< Linear strain
    ST Egl;                    ///< Green-Lagrange strain
    Eigen::VectorX<ST> dEdu;   ///< first derivative of Egl w.r.t displacements
    Eigen::MatrixX<ST> ddEddu; ///< second derivative of Egl w.r.t displacements
  };

  /**
   * \brief Constructor for the Truss class.
   * \param pre The pre fe
   */
  explicit Truss(const Pre& pre)
      : E_{pre.E},
        A_{pre.A},
        density_{pre.density} {}

  /**
   * \brief Gets the displacement for the given Requirement and optional displacement vector.
   *
   * \tparam ST The scalar type for the displacement vector.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \return The local displacement vector of the type Dune::BlockVector<Dune::RealTuple<ST, worldDim>>
   */
  template <typename ST = double>
  auto displacement(const Requirement& par,
                    const std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>& dx = std::nullopt) const {
    const auto& d = par.globalSolution();
    auto disp     = Ikarus::FEHelper::localSolutionBlockVector<Traits>(d, underlying().localView(), dx);
    return disp;
  }

  /**
   * \brief Gets the strain for the given Requirement and optional displacement vector.
   *
   * \tparam ST The scalar type for the strain.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \return A tuple containing Green-Lagrange strain, its first and second derivatives w.r.t displacements and
   * the length of the undeformed and deformed geometry.
   */
  template <class ST = double>
  KinematicVariables<ST> computeStrain(
      const Requirement& par,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>& dx = std::nullopt) const {
    KinematicVariables<ST> kin;
    const auto& localView    = underlying().localView();
    const auto& tree         = localView.tree();
    const auto numberOfNodes = tree.child(0).finiteElement().size();
    auto& ele                = localView.element();
    const auto X1            = Dune::toEigen(ele.geometry().corner(0));
    const auto X2            = Dune::toEigen(ele.geometry().corner(1));
    const auto u             = Dune::viewAsEigenMatrixAsDynFixed(displacement(par, dx)).transpose().eval();
    const auto A1            = X2 - X1;
    const auto ud1           = u.col(1) - u.col(0);

    const Eigen::Vector<ST, worldDim> x1 = X1 + u.col(0);
    const Eigen::Vector<ST, worldDim> x2 = X2 + u.col(1);

    const double Lsquared = A1.squaredNorm();
    const ST lsquared     = (x2 - x1).squaredNorm();

    kin.L  = sqrt(Lsquared);
    kin.l  = sqrt(lsquared);
    kin.A1 = A1;

    // Linear strain
    kin.Elin = (kin.l - kin.L) / kin.L;

    // Green-Lagrange strains
    kin.Egl = 0.5 * (lsquared - Lsquared) / Lsquared;

    // First derivative of Egl w.r.t displacements
    kin.dEdu.setZero(worldDim * numberOfNodes);
    kin.dEdu << -A1 - ud1, A1 + ud1;
    kin.dEdu /= Lsquared;

    // Second derivative of Egl w.r.t displacements
    kin.ddEddu.setIdentity(worldDim * numberOfNodes, worldDim * numberOfNodes);
    kin.ddEddu.template topRightCorner<worldDim, worldDim>()   = -Eigen::Matrix<ST, worldDim, worldDim>::Identity();
    kin.ddEddu.template bottomLeftCorner<worldDim, worldDim>() = -Eigen::Matrix<ST, worldDim, worldDim>::Identity();
    kin.ddEddu /= Lsquared;

    return kin;
  }

  /**
   * \brief Computes the length of the truss and the coordinate transformation matrix T for the undeformed configuration
   *
   * \return std::pair of length, and T
   */
  auto computeLengthAndTransformationMatrix(const Requirement& par) const {
    const auto [L, A1, l, Elin, Egl, dEdu, ddEddu] = computeStrain(par);
    const auto A1normed                            = A1 / L;

    auto T = Eigen::Matrix<typename std::remove_cvref_t<decltype(L)>, 2, worldDim * 2>::Zero().eval();
    T.template block<1, worldDim>(0, 0)        = A1normed;
    T.template block<1, worldDim>(1, worldDim) = A1normed;

    return std::make_pair(L, T);
  }

public:
  /**
   * \brief Calculates a requested result at a specific local position.
   *
   * \param req The Requirement object holding the global solution.
   * \param local Local position vector.
   * \tparam RT The requested result type
   * \return calculated result
   *
   * \tparam RT The type representing the requested result.
   */
  template <template <typename, int, int> class RT>
  requires(supportsResultType<RT>())
  auto calculateAtImpl(const Requirement& req, [[maybe_unused]] const Dune::FieldVector<double, Traits::mydim>& local,
                       Dune::PriorityTag<0>) const {
    using RTWrapper                               = ResultWrapper<RT<double, myDim, myDim>, ResultShape::Vector>;
    const auto [L, _, l, Elin, Egl, dEdu, ddEddu] = computeStrain(req);
    if constexpr (isSameResultType<RT, ResultTypes::cauchyAxialForce>) {
      auto N = Eigen::Vector<double, 1>{E_ * A_ * Egl * l / L}; // Axial force in deformed configuration
      return RTWrapper{N};
    }
    if constexpr (isSameResultType<RT, ResultTypes::PK2AxialForce>) {
      auto N = Eigen::Vector<double, 1>{E_ * A_ * Egl}; // Axial force in undeformed configuration
      return RTWrapper{N};
    }
    if constexpr (isSameResultType<RT, ResultTypes::linearAxialForce>) {
      auto N = Eigen::Vector<double, 1>{E_ * A_ * Elin};
      return RTWrapper{N};
    }
  }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const FE&>(*this); }
  auto& underlying() { return static_cast<FE&>(*this); }

  double E_;
  double A_;
  double density_;

protected:
  template <typename ScalarType>
  void calculateMatrixImpl(
      const Requirement& par, const MatrixAffordance& affordance, typename Traits::template MatrixType<> K,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    if (affordance == MatrixAffordance::stiffness) {
      const auto [L, _, l, Elin, Egl, dEdu, ddEddu] = computeStrain(par, dx);
      K += E_ * A_ * L * (dEdu * dEdu.transpose() + ddEddu * Egl);
    } else if (affordance == MatrixAffordance::linearMass) {
      Eigen::Matrix<ScalarType, 2, 2> mLoc{
          {2, 1},
          {1, 2}
      };
      auto [l, T] = computeLengthAndTransformationMatrix(par);
      mLoc *= (density_ * A_ * l) / 6.0;
      K += T.transpose() * mLoc * T;
    } else
      DUNE_THROW(Dune::NotImplemented, "MatrixAffordance not implemented: " + toString(affordance));
  }

  template <typename ScalarType>
  auto calculateScalarImpl(const Requirement& par, ScalarAffordance affordance,
                           const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx =
                               std::nullopt) const -> ScalarType {
    const auto [L, l, _, Elin, Egl, dEdu, ddEddu] = computeStrain(par, dx);
    return 0.5 * E_ * A_ * L * Egl * Egl;
  }

  template <typename ScalarType>
  void calculateVectorImpl(
      const Requirement& par, VectorAffordance affordance, typename Traits::template VectorType<ScalarType> force,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    const auto [L, l, _, Elin, Egl, dEdu, ddEddu] = computeStrain(par, dx);
    force += E_ * A_ * Egl * L * dEdu;
  }
};

/**
 * \brief A helper function to create a truss pre finite element.
 * \param E Young's modulus of the truss member
 * \param A Cross section area of the truss member
 * \param density Density of truss member (defaults to 1.0)
 * \return A truss pre finite element.
 */
inline auto truss(const double E, const double A, double density = 1.0) {
  TrussPre pre(E, A, density);

  return pre;
}
} // namespace Ikarus
