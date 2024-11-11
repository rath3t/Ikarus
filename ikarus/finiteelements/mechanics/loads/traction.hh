// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/fufem/boundarypatch.hh>
#include <dune/localfefunctions/derivativetransformators.hh>
#include <dune/localfefunctions/meta.hh>

#include <ikarus/finiteelements/ferequirements.hh>

namespace Ikarus {

template <typename PreFE, typename FE>
class Traction;

/**
 * \brief A PreFE struct for Neumann boundary load skill.
 * \tparam GV Type of the grid view.
 */
template <typename GV>
struct NeumannBoundaryLoadPre
{
  using GridView                = GV;
  static constexpr int worldDim = GridView::dimensionworld;
  const BoundaryPatch<GridView>* neumannBoundary;

  std::function<Eigen::Vector<double, worldDim>(const Dune::FieldVector<double, worldDim>&, const double&)> load;
  using BoundaryPatchType = BoundaryPatch<GridView>;

  template <typename PreFE, typename FE>
  using Skill = Traction<PreFE, FE>;
};

/**
 * \brief Traction class represents distributed traction load that can be applied.
 *
 * \ingroup mechanics
 *
 * \tparam PreFE The type of the  pre finite element.
 * \tparam FE The type of the finite element.
 */
template <typename PreFE, typename FE>
class Traction
{
public:
  using Traits = PreFE::Traits;
  using Requirement =
      FERequirementsFactory<FESolutions::displacement, FEParameter::loadfactor, Traits::useEigenRef>::type;
  using LocalView               = typename Traits::LocalView;
  using GridView                = typename Traits::GridView;
  static constexpr int myDim    = Traits::mydim;
  static constexpr int worldDim = Traits::worlddim;
  using Pre                     = NeumannBoundaryLoadPre<GridView>;

  /**
   * \brief Constructor for the traction class.
   *
   * \param pre The pre finite element
   */
  explicit Traction(const Pre& pre)
      : neumannBoundaryLoad_{pre.load},
        neumannBoundary_{pre.neumannBoundary} {}

protected:
  template <template <typename, int, int> class RT>
  requires Dune::AlwaysFalse<RT<double, 1, 1>>::value
  auto calculateAtImpl(const Requirement& req, const Dune::FieldVector<double, Traits::mydim>& local,
                       Dune::PriorityTag<0>) const {}

  template <typename ST>
  auto calculateScalarImpl(
      const Requirement& par, ScalarAffordance affordance,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>& dx = std::nullopt) const -> ST {
    if (not neumannBoundary_ and not neumannBoundaryLoad_)
      return 0.0;
    ST energy            = 0.0;
    const auto uFunction = underlying().displacementFunction(par, dx);
    const auto& lambda   = par.parameter();
    auto& element        = underlying().localView().element();

    for (auto&& intersection : intersections(neumannBoundary_->gridView(), element)) {
      if (not neumannBoundary_->contains(intersection))
        continue;

      const auto& quadLine = Dune::QuadratureRules<double, myDim - 1>::rule(intersection.type(), underlying().order());

      for (const auto& curQuad : quadLine) {
        // Local position of the quadrature point
        const Dune::FieldVector<double, myDim>& quadPos = intersection.geometryInInside().global(curQuad.position());

        const double integrationElement = intersection.geometry().integrationElement(curQuad.position());

        // The value of the local function
        const auto uVal = uFunction.evaluate(quadPos);

        // Value of the Neumann data at the current position
        auto neumannValue = neumannBoundaryLoad_(intersection.geometry().global(curQuad.position()), lambda);

        energy -= neumannValue.dot(uVal) * curQuad.weight() * integrationElement;
      }
    }
    return energy;
  }

  template <typename ST>
  void calculateVectorImpl(
      const Requirement& par, VectorAffordance affordance, typename Traits::template VectorType<ST> force,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>& dx = std::nullopt) const {
    if (not neumannBoundary_ and not neumannBoundaryLoad_)
      return;
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto uFunction = underlying().displacementFunction(par, dx);
    const auto& lambda   = par.parameter();
    auto& element        = underlying().localView().element();

    for (auto&& intersection : intersections(neumannBoundary_->gridView(), element)) {
      if (not neumannBoundary_->contains(intersection))
        continue;

      /// Integration rule along the boundary
      const auto& quadLine = Dune::QuadratureRules<double, myDim - 1>::rule(intersection.type(), underlying().order());

      for (const auto& curQuad : quadLine) {
        const Dune::FieldVector<double, myDim>& quadPos = intersection.geometryInInside().global(curQuad.position());
        const double intElement = intersection.geometry().integrationElement(curQuad.position());

        /// The value of the local function wrt the i-th coeff
        for (size_t i = 0; i < underlying().numberOfNodes(); ++i) {
          const auto udCi = uFunction.evaluateDerivative(quadPos, wrt(coeff(i)));

          /// Value of the Neumann data at the current position
          auto neumannValue = neumannBoundaryLoad_(intersection.geometry().global(curQuad.position()), lambda);
          force.template segment<worldDim>(worldDim * i) -= udCi * neumannValue * curQuad.weight() * intElement;
        }
      }
    }
  }

  template <typename ST>
  void calculateMatrixImpl(
      const Requirement&, MatrixAffordance, typename Traits::template MatrixType<>,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>& = std::nullopt) const {}

  void updateStateImpl(const Requirement& /* par */,
                       typename Traits::template VectorTypeConst<> /* correction */) const {}

private:
  std::function<Eigen::Vector<double, worldDim>(const Dune::FieldVector<double, worldDim>&, const double&)>
      neumannBoundaryLoad_;
  const BoundaryPatch<GridView>* neumannBoundary_;

  //> CRTP
  const auto& underlying() const { return static_cast<const FE&>(*this); }
  auto& underlying() { return static_cast<FE&>(*this); }
};

/**
 * \brief A helper function to create a Neumann boundary load skill.
 * \tparam GV Type of the grid view.
 * \tparam F Type of the Neumann boundary load functor.
 * \param patch The patch where Neumann boundary load is applied.
 * \param load The neumann boundary load functor.
 * \return A Neumann boundary load skill.
 */
template <typename GV, typename F>
auto neumannBoundaryLoad(const BoundaryPatch<GV>* patch, F&& load) {
  NeumannBoundaryLoadPre<GV> pre(patch, std::forward<F>(load));

  return pre;
}

} // namespace Ikarus
