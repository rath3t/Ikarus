// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file interface.hh
 * \brief Implementation of the Hyperelastic material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/concepts.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/volumetric/volumetricfunctions.hh>
#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Implementation of a general Hyperelastic Material material model.
 * \details \f$\Psi(\BC) = \hat{\Psi}(\la_1, \la_2, \la_3) + U(J)\f$ with \f$\hat{\Psi}\f$ being the
 *deviatoric part of the strain energy function and \f$ U(J) \f$ being the volumetric part. After calling the underlying
 *deviatoric and volumetric function, the transformation to cartesian coordinate system is implemented in this
 *interface.
 *
 * \ingroup materials
 */
template <typename DEV, typename VOL = NoVolumetricPart>
requires(std::same_as<typename DEV::ScalarType, typename VOL::ScalarType>)
struct Hyperelastic : public Material<Hyperelastic<DEV, VOL>>
{
  // Checking concepts here results in better compiler error messages because of CRTP
  static_assert(Concepts::DeviatoricPart<DEV>);
  static_assert(Concepts::VolumetricPart<VOL>);

  using ScalarType                        = typename DEV::ScalarType;
  static constexpr bool hasVolumetricPart = not std::same_as<VOL, NoVolumetricPart>;

  static constexpr int dim = 3;
  using StrainMatrix       = Eigen::Matrix<ScalarType, dim, dim>;
  using StressMatrix       = StrainMatrix;
  using MaterialTensor     = Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>>;

  using MaterialParametersDEV = typename DEV::MaterialParameters;
  using MaterialParametersVOL = typename VOL::MaterialParameter;
  using MaterialParameters =
      std::conditional_t<hasVolumetricPart, std::pair<MaterialParametersDEV, MaterialParametersVOL>,
                         MaterialParametersDEV>;

  static constexpr auto strainTag              = StrainTags::rightCauchyGreenTensor;
  static constexpr auto stressTag              = StressTags::PK2;
  static constexpr auto tangentModuliTag       = TangentModuliTags::Material;
  static constexpr bool energyAcceptsVoigt     = true;
  static constexpr bool stressToVoigt          = false;
  static constexpr bool stressAcceptsVoigt     = true;
  static constexpr bool moduliToVoigt          = false;
  static constexpr bool moduliAcceptsVoigt     = true;
  static constexpr double derivativeFactorImpl = 2;

  [[nodiscard]] constexpr static std::string nameImpl() noexcept {
    if constexpr (hasVolumetricPart)
      return "Hyperelastic (" + DEV::name() + ", " + VOL::name() + ")";
    else
      return "Hyperelastic (" + DEV::name() + ")";
  }

  explicit Hyperelastic(const DEV& dev)
  requires(not hasVolumetricPart)
      : dev_{dev},
        vol_(NoVolumetricPart{MaterialParametersVOL{}, {}}) {}

  Hyperelastic(const DEV& dev, const VOL& vol)
      : dev_(dev),
        vol_(vol) {}

  /** \brief Returns the deviatoric function. */
  const DEV& deviatoricFunction() const { return dev_; }

  /** \brief Returns the volumetric function. */
  const VOL& volumetricFunction() const { return vol_; }

  /**
   * \brief Returns the material parameters stored in the deviatoric part of the material.
   */
  const MaterialParameters materialParametersImpl() const {
    if constexpr (hasVolumetricPart)
      return {dev_.materialParameters(), vol_.materialParameter()};
    else
      return dev_.materialParameters();
  }

  /**
   * \brief Computes the total stored energy in the Hyperelastic material model.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);

    const auto& lambdas = principalStretches(C, Eigen::EigenvaluesOnly).first;
    auto J              = detF(lambdas);

    return dev_.storedEnergy(lambdas) + vol_.storedEnergy(J);
  }

  /**
   * \brief Computes the stresses in the Hyperelastic material model.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return StressMatrix The stresses.
   */
  template <bool voigt, typename Derived>
  StressMatrix stressesImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      const auto& [lambdas, N] = principalStretches(C);
      auto J                   = detF(lambdas);

      const auto& Sdev = transformDeviatoricStresses(dev_.stresses(lambdas), N);
      const auto& Svol = transformVolumetricStresses(vol_.firstDerivative(J), C, J);

      return Sdev + Svol;
    } else
      static_assert(voigt == false, "Hyperelastic does not support returning stresses in Voigt notation");
  }

  /**
   * \brief Computes the tangent moduli in the Hyperelastic material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> The tangent moduli.
   */
  template <bool voigt, typename Derived>
  MaterialTensor tangentModuliImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      const auto& [lambdas, N] = principalStretches(C);
      auto J                   = detF(lambdas);

      const auto& moduliDev = transformDeviatoricTangentModuli(dev_.tangentModuli(lambdas), N);
      const auto& moduliVol = transformVolumetricTangentModuli(vol_.firstDerivative(J), vol_.secondDerivative(J), C, J);

      return moduliDev + moduliVol;
    } else
      static_assert(voigt == false, "Hyperelastic does not support returning tangent moduli in Voigt notation");
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return Hyperelastic<ScalarTypeOther> The rebound Hyperelastic material.
   */
  template <typename STO>
  auto rebind() const {
    auto reboundDEV = dev_.template rebind<STO>();
    auto reboundVOL = vol_.template rebind<STO>();
    return Hyperelastic<decltype(reboundDEV), decltype(reboundVOL)>(reboundDEV, reboundVOL);
  }

private:
  DEV dev_;
  VOL vol_;

  inline auto dimensionRange() const { return Dune::range(dim); }

  StressMatrix transformDeviatoricStresses(const typename DEV::StressMatrix& principalStress,
                                           const Eigen::Matrix<ScalarType, 3, 3>& N) const {
    auto S = StressMatrix::Zero().eval();
    for (auto i : Dune::range(3))
      S += principalStress[i] * dyadic<ScalarType, 3, false>(N.col(i).eval(), N.col(i).eval());

    return S;
  }

  StressMatrix transformVolumetricStresses(const ScalarType& Uprime, const auto& C, ScalarType J) const {
    return J * Uprime * C.inverse();
  }

  MaterialTensor transformDeviatoricTangentModuli(const typename DEV::MaterialTensor& L,
                                                  const Eigen::Matrix<ScalarType, 3, 3>& N) const {
    MaterialTensor moduli{};
    moduli.setZero();

    for (const auto i : dimensionRange())
      for (const auto k : dimensionRange()) {
        // First term: L[i, i, k, k] * ((N[i] ⊗ N[i]) ⊗ (N[k] ⊗ N[k]))
        auto NiNi = dyadic(N.col(i).eval(), N.col(i).eval());
        auto NkNk = dyadic(N.col(k).eval(), N.col(k).eval());

        moduli += L(i, i, k, k) * dyadic(NiNi, NkNk);

        // Second term (only if i != k): L[i, k, i, k] * (N[i] ⊗ N[k] ⊗ (N[i] ⊗ N[k] + N[k] ⊗ N[i]))
        if (i != k) {
          auto NiNk = dyadic(N.col(i).eval(), N.col(k).eval());
          auto NkNi = dyadic(N.col(k).eval(), N.col(i).eval());

          moduli += L(i, k, i, k) * dyadic(NiNk, NiNk + NkNi);
        }
      }

    return moduli;
  }

  MaterialTensor transformVolumetricTangentModuli(const ScalarType& Uprime, const ScalarType& Uprimeprime,
                                                  const auto& C, ScalarType J) const {
    const auto invC  = C.inverse().eval();
    const auto CTinv = tensorView(invC, std::array<Eigen::Index, 2>({3, 3}));

    MaterialTensor moduli = (J * ((Uprime + J * Uprimeprime) * dyadic(CTinv, CTinv) -
                                  (2 * Uprime * symTwoSlots(fourthOrderIKJL(invC, invC), {2, 3}))))
                                .eval();

    return moduli;
  }

  template <typename Derived>
  auto principalStretches(const Eigen::MatrixBase<Derived>& Craw, int options = Eigen::ComputeEigenvectors) const {
    StrainMatrix C = Impl::maybeFromVoigt(Craw);
    return Impl::principalStretches<ScalarType>(C, options);
  }

  auto detF(const typename DEV::PrincipalStretches& lambda) const -> typename VOL::JType {
    if constexpr (hasVolumetricPart) {
      const auto detC = Impl::determinantFromPrincipalValues<ScalarType>(lambda);
      Impl::checkPositiveDet(detC);
      return detC;
    }
    return 0.0;
  }
};

} // namespace Ikarus::Materials
