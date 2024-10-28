// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file BlatzKo.hh
 * \brief Implementation of the BlatzKo material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

/**
 * \brief Implementation of the Neo-Hookean material model.
* \ingroup materials
*  The energy is computed as
*  \f[ \psi(\BC) = \frac{\mu}{2} (\tr \BC-3- 2 \log \sqrt{\det \BC}) + \frac{\lambda}{2} (\log \sqrt{\det \BC})^2 ,\f]
* where \f$ \BC \f$ denotes the right Cauchy-Green strain tensor.
*
*  The second Piola-Kirchhoff stresses are computed as
*      \f[ \BS(\BC) =\fracpt{\psi(\BC)}{\BC} = \mu (\BI-\BC^{-1}) + \lambda \log \sqrt{\det \BC}  \BC^{-1},\f]
*
* and the material tangent moduli are computed as
*      \f[ \BBC(\BC) =\fracpt{^2\psi(\BC)}{\BC^2} =  \lambda \BC^{-1} \otimes  \BC^{-1} + 2 (\mu- \lambda \log
\sqrt{\det \BC} ) \CI,\f]
*      where \f$ \CI_{IJKL} =  \frac{1}{2}({(\BC^{-1})}^{IK}{(\BC^{-1})}^{JL}+{(\BC^{-1})}^{IL} {(\BC^{-1})}^{JK}).\f$
*
*  \remark See \cite bonet2008nonlinear, Section 6.4.3 for a discussion of this material
 * \tparam ST The scalar type for the strains and stresses,....
 */
template <typename ST>
struct BlatzKoT : public Material<BlatzKoT<ST>>
{
  using ScalarType                    = ST;
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

  [[nodiscard]] constexpr static std::string nameImpl() noexcept { return "BlatzKo"; }

  /**
   * \brief Constructor for BlatzKoT.
   * \param mpt The Lame's parameters (first parameter and shear modulus).
   */
  explicit BlatzKoT(const MaterialParameters& mpt)
      : materialParameter_{mpt} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameter_; }

  template <typename Derived>
  auto principalStretches(const Eigen::MatrixBase<Derived>& C) const {
    Eigen::SelfAdjointEigenSolver<Derived> eigensolver(C);
    auto& eigenvalues  = eigensolver.eigenvalues();
    auto& eigenvectors = eigensolver.eigenvectors();

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
      auto [lambdas, N] = principalStretches(C);

      return materialParameter_.mu / 2 *
             (1 / std::pow(lambdas[0], 2) + 1 / std::pow(lambdas[1], 2) + 1 / std::pow(lambdas[2], 2) +
              2 * lambdas[0] * lambdas[1] * lambdas[2] - 5);

    } else
      static_assert(!Concepts::EigenVector<Derived>,
                    "BlatzKo energy can only be called with a matrix and not a vector in Voigt notation");
  }

  auto principalStresseses(const auto& lambdas) const {
    // principal PK1 stress
    auto P1 = materialParameter_.mu * (-2 / std::pow(lambdas[0], 3) + 2 * lambdas[1] * lambdas[2]) / 2;
    auto P2 = materialParameter_.mu * (-2 / std::pow(lambdas[1], 3) + 2 * lambdas[0] * lambdas[2]) / 2;
    auto P3 = materialParameter_.mu * (-2 / std::pow(lambdas[2], 3) + 2 * lambdas[0] * lambdas[1]) / 2;

    // principal PK2 stress
    auto S1 = 1 / lambdas[0] * P1;
    auto S2 = 1 / lambdas[1] * P2;
    auto S3 = 1 / lambdas[2] * P3;

    return Eigen::Vector<ScalarType, 3>{S1, S2, S3};
  }

  template <typename ST_, int size>
  auto dyadicProduct(const Eigen::Vector<ST_, size>& a, const Eigen::Vector<ST_, size>& b) const {
    return (a * b.transpose()).eval();
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
        auto [lambdas, N] = principalStretches(C);

        auto principalStress = principalStresseses(lambdas);

        // Transformation from principal coordinates to cartesian coordinates
        auto S = Eigen::Matrix3<ScalarType>::Zero().eval();
        for (auto i : Dune::range(3))
          S += principalStress[i] * dyadicProduct(N.col(i).eval(), N.col(i).eval());

        return S;
      } else
        static_assert(!Concepts::EigenVector<Derived>,
                      "BlatzKo can only be called with a matrix and not a vector in Voigt notation");
    } else
      static_assert(voigt == false, "BlatzKo does not support returning stresses in Voigt notation");
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
      const auto invC = C.inverse().eval();
      const auto detC = C.determinant();
      checkPositiveDetC(detC);
      const auto logdetF = log(sqrt(detC));
      const auto CTinv   = tensorView(invC, std::array<Eigen::Index, 2>({3, 3}));
      static_assert(Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3>>::NumIndices == 2);
      Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> moduli =
          (materialParameter_.lambda * dyadic(CTinv, CTinv) +
           2 * (materialParameter_.mu - materialParameter_.lambda * logdetF) *
               symTwoSlots(fourthOrderIKJL(invC, invC), {2, 3}))
              .eval();
      return moduli;
    } else
      static_assert(voigt == false, "BlatzKo does not support returning tangent moduli in Voigt notation");
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return BlatzKoT<ScalarTypeOther> The rebound BlatzKo material.
   */
  template <typename STO>
  auto rebind() const {
    return BlatzKoT<STO>(materialParameter_);
  }

private:
  MaterialParameters materialParameter_;

  // TODO reimplemet
  void checkPositiveDetC(ScalarType detC) const {
    if (Dune::FloatCmp::le(static_cast<double>(detC), 0.0, 1e-10))
      DUNE_THROW(Dune::InvalidStateException,
                 "Determinant of right Cauchy Green tensor C must be greater than zero. detC = " +
                     std::to_string(static_cast<double>(detC)));
  }
};

/**
 * \brief Alias for BlatzKoT with double as the default scalar type.
 */
using BlatzKo = BlatzKoT<double>;

} // namespace Ikarus
