// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/localfefunctions/derivativetransformators.hh>
#include <dune/localfefunctions/meta.hh>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {

template <typename PreFE, typename FE>
class VolumeLoad;

/**
 * \brief A PreFE struct for volume load skill.
 * \tparam wd The world dimension.
 */
template <int wd>
struct VolumeLoadPre
{
  static constexpr int worldDim = wd;
  std::function<Eigen::Vector<double, worldDim>(const Dune::FieldVector<double, worldDim>&, const double&)> volumeLoad;

  template <typename PreFE, typename FE>
  using Skill = VolumeLoad<PreFE, FE>;
};

// Deduction guide
template <class F>
VolumeLoadPre(F f) -> VolumeLoadPre<traits::FunctionTraits<F>::return_type::RowsAtCompileTime>;

/**
 * \brief VolumeLoad class represents distributed volume load that can be applied.
 *
 * \ingroup mechanics
 *
 * \tparam PreFE The type of the  pre finite element.
 * \tparam FE The type of the finite element.
 */
template <typename PreFE, typename FE>
class VolumeLoad
{
public:
  using Traits = PreFE::Traits;
  using Requirement =
      FERequirementsFactory<FESolutions::displacement, FEParameter::loadfactor, Traits::useEigenRef>::type;
  static constexpr int worldDim = Traits::worlddim;
  using Pre                     = VolumeLoadPre<worldDim>;

  /**
   * \brief Constructor for the Loads class.
   *
   * \param pre Volume load pre object.
   */
  explicit VolumeLoad(const Pre& pre = {})
      : volumeLoad_{pre.volumeLoad} {}

protected:
  template <template <typename, int, int> class RT>
  requires Dune::AlwaysFalse<RT<double, 1, 1>>::value
  auto calculateAtImpl(const Requirement& req, const Dune::FieldVector<double, Traits::mydim>& local,
                       Dune::PriorityTag<0>) const {}

  template <typename ST>
  auto calculateScalarImpl(
      const Requirement& par, ScalarAffordance affordance,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>& dx = std::nullopt) const -> ST {
    if (not volumeLoad_)
      return 0.0;
    ST energy            = 0.0;
    const auto uFunction = underlying().displacementFunction(par, dx);
    const auto geo       = underlying().geometry();
    const auto& lambda   = par.parameter();

    for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
      const auto uVal                      = uFunction.evaluate(gpIndex);
      Eigen::Vector<double, worldDim> fext = volumeLoad_(geo.global(gp.position()), lambda);
      energy -= uVal.dot(fext) * geo.integrationElement(gp.position()) * gp.weight();
    }
    return energy;
  }

  template <typename ST>
  void calculateVectorImpl(
      const Requirement& par, VectorAffordance affordance, typename Traits::template VectorType<ST> force,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>& dx = std::nullopt) const {
    if (not volumeLoad_)
      return;
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto uFunction = underlying().displacementFunction(par, dx);
    const auto geo       = underlying().geometry();
    const auto& lambda   = par.parameter();

    for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
      const Eigen::Vector<double, worldDim> fext = volumeLoad_(geo.global(gp.position()), lambda);
      const double intElement                    = geo.integrationElement(gp.position()) * gp.weight();
      for (size_t i = 0; i < underlying().numberOfNodes(); ++i) {
        const auto udCi = uFunction.evaluateDerivative(gpIndex, wrt(coeff(i)));
        force.template segment<worldDim>(worldDim * i) -= udCi * fext * intElement;
      }
    }
  }

  template <typename ST>
  void calculateMatrixImpl(
      const Requirement& par, MatrixAffordance, typename Traits::template MatrixType<> K,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>& dx = std::nullopt) const {}

  void updateStateImpl(const Requirement& /* par */,
                       typename Traits::template VectorTypeConst<> /* correction */) const {}

private:
  std::function<Eigen::Vector<double, worldDim>(const Dune::FieldVector<double, worldDim>&, const double&)> volumeLoad_;
  //> CRTP
  const auto& underlying() const { return static_cast<const FE&>(*this); }
  auto& underlying() { return static_cast<FE&>(*this); }
};

/**
 * \brief A helper function to create a volume load skill.
 * \tparam worldDim The world dimension.
 * \param f A function representing the volume load.
 * \return A volume load skill.
 */
template <int worldDim>
auto volumeLoad(const std::function<Eigen::Vector<double, worldDim>(const Dune::FieldVector<double, worldDim>&,
                                                                    const double&)>& f) {
  VolumeLoadPre pre(f);
  return pre;
}

/**
 * \brief A helper function to create a volume load skill.
 * \param f A function representing the volume load.
 * \return A volume load skill.
 */
template <typename F>
auto volumeLoad(F&& f) {
  VolumeLoadPre pre(f);
  return pre;
}

} // namespace Ikarus
