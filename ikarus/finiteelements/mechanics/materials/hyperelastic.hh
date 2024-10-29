// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file Hyperelastic.hh
 * \brief Implementation of the Hyperelastic material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

/**
 * \brief Implementation of a general Hyperelastic Material material model.
 * \ingroup materials
 */
template <typename DEV>
struct Hyperelastic : public Material<Hyperelastic<DEV>>
{
  using ScalarType     = typename DEV::ScalarType;
  using DeviatoricPart = DEV;

  static constexpr int worldDimension = 3;
  using StrainMatrix                  = Eigen::Matrix<ScalarType, worldDimension, worldDimension>;
  using StressMatrix                  = StrainMatrix;
  using MaterialParameters            = LamesFirstParameterAndShearModulus;

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

  /**
   * \brief Constructor for Hyperelastic.
   * \param mpt The Lame's parameters (first parameter and shear modulus).
   */
  explicit Hyperelastic(const DEV& dev)
      : dev_{dev} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return dev_.materialParameter_; }

  // TODO: Can we use SelfAdjointSolver?
  template <typename Derived>
  auto principalStretches(const Eigen::MatrixBase<Derived>& C, int options = Eigen::ComputeEigenvectors) const {
    Eigen::SelfAdjointEigenSolver<Derived> eigensolver(C, options);
    auto& eigenvalues  = eigensolver.eigenvalues();
    auto& eigenvectors = options == Eigen::ComputeEigenvectors ? eigensolver.eigenvectors() : Derived::Zero();

    auto principalStretches = eigenvalues.array().sqrt().eval();
    return std::make_pair(principalStretches, eigenvectors);
  }

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
      return dev_.storedEnergyImpl(lambdas);

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
  auto stressesImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      if constexpr (!Concepts::EigenVector<Derived>) {
        auto [lambdas, N]    = principalStretches(C);
        auto principalStress = dev_.stressesImpl(lambdas);

        // Transformation from principal coordinates to cartesian coordinates
        auto S = Eigen::Matrix3<ScalarType>::Zero().eval();
        for (auto i : Dune::range(3))
          S += principalStress[i] * dyadic(N.col(i).eval(), N.col(i).eval());

        return S;
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
  auto tangentModuliImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      auto [lambdas, N] = principalStretches(C);

      auto L = dev_.tangentModuliImpl(lambdas);

      Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> moduli{};
      moduli.setZero();

      for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
          // First term: L[i, i, k, k] * ((N[i] ⊗ N[i]) ⊗ (N[k] ⊗ N[k]))
          auto NiNi = tensorView(dyadic(N.col(i).eval(), N.col(i).eval()), std::array<Eigen::Index, 2>({3, 3}));
          auto NkNk = tensorView(dyadic(N.col(k).eval(), N.col(k).eval()), std::array<Eigen::Index, 2>({3, 3}));

          moduli += L(i, i, k, k) * dyadic(NiNi, NkNk);

          // Second term (only if i != k): L[i, k, i, k] * (N[i] ⊗ N[k] ⊗ (N[i] ⊗ N[k] + N[k] ⊗ N[i]))
          if (i != k) {
            auto NiNk = tensorView(dyadic(N.col(i).eval(), N.col(k).eval()), std::array<Eigen::Index, 2>({3, 3}));
            auto NkNi = tensorView(dyadic(N.col(k).eval(), N.col(i).eval()), std::array<Eigen::Index, 2>({3, 3}));

            moduli += L(i, k, i, k) * dyadic(NiNk, NiNk + NkNi);
          }
        }
      }
      return moduli;
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
    return Hyperelastic<decltype(reboundDEV)>(reboundDEV);
  }

private:
  DEV dev_;
};

} // namespace Ikarus
