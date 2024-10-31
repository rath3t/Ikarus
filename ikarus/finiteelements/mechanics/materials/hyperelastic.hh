// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file Hyperelastic.hh
 * \brief Implementation of the Hyperelastic material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/volumetric.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

/**
 * \brief Implementation of a general Hyperelastic Material material model.
 * \ingroup materials
 */
template <typename DEV, typename VOL = NoVolumetricPart>
requires(std::same_as<typename DEV::ScalarType, typename VOL::ScalarType>)
struct Hyperelastic : public Material<Hyperelastic<DEV, VOL>>
{
  using ScalarType                        = typename DEV::ScalarType;
  static constexpr bool hasVolumetricPart = not std::same_as<VOL, NoVolumetricPart>;

  static constexpr int worldDimension = 3;
  using StrainMatrix                  = Eigen::Matrix<ScalarType, worldDimension, worldDimension>;
  using StressMatrix                  = StrainMatrix;
  using MaterialTensor                = Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>>;

  using MaterialParametersDEV = typename DEV::MaterialParameters;
  using MaterialParameters =
      std::conditional_t<hasVolumetricPart, std::pair<MaterialParametersDEV, BulkModulus>, MaterialParametersDEV>;

  static constexpr auto strainTag              = StrainTags::rightCauchyGreenTensor;
  static constexpr auto stressTag              = StressTags::PK2;
  static constexpr auto tangentModuliTag       = TangentModuliTags::Material;
  static constexpr bool energyAcceptsVoigt     = false;
  static constexpr bool stressToVoigt          = false;
  static constexpr bool stressAcceptsVoigt     = false;
  static constexpr bool moduliToVoigt          = false;
  static constexpr bool moduliAcceptsVoigt     = false;
  static constexpr double derivativeFactorImpl = 2;

  [[nodiscard]] constexpr static std::string nameImpl() noexcept { return "Hyperelastic"; }

  explicit Hyperelastic(const DEV& dev)
  requires(not hasVolumetricPart)
      : dev_{dev},
        vol_(NoVolumetricPart{}) {}

  Hyperelastic(const DEV& dev, const VOL& vol)
      : dev_(dev),
        vol_(vol) {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return dev_.materialParameter_; }

  /**
   * \brief Computes the stored energy in the Neo-Hookean material model.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!Concepts::EigenVector<Derived>) {
      auto lambdas = principalStretches(C, Eigen::EigenvaluesOnly).first;
      auto J       = detF(C);

      return dev_.storedEnergyImpl(lambdas) + vol_.storedEnergy(J);

    } else
      static_assert(!Concepts::EigenVector<Derived>,
                    "Hyperelastic energy can only be called with a matrix and not a vector in Voigt notation");
  }

  /**
   * \brief Computes the stresses in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return StressMatrix The stresses.
   */
  template <bool voigt, typename Derived>
  StressMatrix stressesImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      if constexpr (!Concepts::EigenVector<Derived>) {
        auto [lambdas, N] = principalStretches(C);
        auto J            = detF(C);

        auto Sdev = transformDeviatoricStresses(dev_.stressesImpl(lambdas), N);
        auto Svol = transformVolumetricStresses(vol_.firstDerivative(J), C, J);

        return Sdev + Svol;
      } else
        static_assert(!Concepts::EigenVector<Derived>,
                      "Hyperelastic can only be called with a matrix and not a vector in Voigt notation");
    } else
      static_assert(voigt == false, "Hyperelastic does not support returning stresses in Voigt notation");
  }

  /**
   * \brief Computes the tangent moduli in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> The tangent moduli.
   */
  template <bool voigt, typename Derived>
  MaterialTensor tangentModuliImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      auto [lambdas, N] = principalStretches(C);
      auto J            = detF(C);

      auto moduliDev = transformDeviatoricTangentModuli(dev_.tangentModuliImpl(lambdas), N);
      auto moduliVol = transformVolumetricTangentModuli(vol_.firstDerivative(J), vol_.secondDerivative(J), C, J);

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
    auto reboundVOL = dev_.template rebind<VOL>();
    return Hyperelastic<decltype(reboundDEV), decltype(reboundVOL)>(reboundDEV);
  }

private:
  DEV dev_;
  VOL vol_;

  StressMatrix transformDeviatoricStresses(const typename DEV::StressMatrix& principalStress, const auto& N) const {
    auto S = StressMatrix::Zero().eval();
    for (auto i : Dune::range(3))
      S += principalStress[i] * dyadic<ScalarType, 3, false>(N.col(i).eval(), N.col(i).eval());

    return S;
  }

  StressMatrix transformVolumetricStresses(const ScalarType& Uprime, const auto& C, ScalarType J) const {
    return J * Uprime * C.inverse();
  }

  MaterialTensor transformDeviatoricTangentModuli(const typename DEV::MaterialTensor& L, const auto& N) const {
    MaterialTensor moduli{};
    moduli.setZero();

    for (int i = 0; i < worldDimension; ++i) {
      for (int k = 0; k < worldDimension; ++k) {
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
    }
    return moduli;
  }

  MaterialTensor transformVolumetricTangentModuli(const ScalarType& Uprime, const ScalarType& Uprimeprime,
                                                  const auto& C, ScalarType J) const {
    const auto invC    = C.inverse().eval();
    const auto CTinv   = tensorView(invC, std::array<Eigen::Index, 2>({3, 3}));
    const auto CinvDya = dyadic(CTinv, CTinv);
    const auto CinvT23 = symTwoSlots(fourthOrderIKJL(invC, invC), {2, 3});

    Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> moduli =
        (J * ((Uprime + J * Uprimeprime) * CinvDya - 2 * Uprime * CinvT23)).eval();

    return moduli;
  }

  // TODO: Can we use SelfAdjointSolver?
  template <typename Derived>
  auto principalStretches(const Eigen::MatrixBase<Derived>& C, int options = Eigen::ComputeEigenvectors) const {
    Eigen::SelfAdjointEigenSolver<Derived> eigensolver(C, options);
    // Eigen::EigenSolver<Derived> eigensolver(C, options);
    auto& eigenvalues  = eigensolver.eigenvalues();
    auto& eigenvectors = options == Eigen::ComputeEigenvectors ? eigensolver.eigenvectors() : Derived::Zero();

    auto principalStretches = eigenvalues.array().sqrt().eval();
    return std::make_pair(principalStretches, eigenvectors);
  }

  template <typename Derived>
  auto detF(const Eigen::MatrixBase<Derived>& C) const -> typename VOL::JType {
    if constexpr (hasVolumetricPart) {
      const auto detC = C.determinant();
      const auto J     = std::sqrt(detC);
      checkPositiveDetC(J);

      return J;
    }
    return 0.0;
  }

  void checkPositiveDetC(ScalarType detC) const {
    if (Dune::FloatCmp::le(static_cast<double>(detC), 0.0, 1e-10))
      DUNE_THROW(Dune::InvalidStateException,
                 "Determinant of right Cauchy Green tensor C must be greater than zero. detC = " +
                     std::to_string(static_cast<double>(detC)));
  }
};

} // namespace Ikarus
