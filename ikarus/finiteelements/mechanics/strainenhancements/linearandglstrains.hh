// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file linearandglstrains.hh
 * \brief Definition of the types of EAS formulations for 2D and 3D linear solid elements, where linear and
 * Green-Lagrange strains are enhanced.
 *
 * \ingroup  eas
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::EAS {

namespace Impl {
  template <typename GEO>
  auto transformationMatrixAtCenterWithDetJ(const GEO& geometry) {
    const auto& referenceElement = Dune::ReferenceElements<double, GEO::mydimension>::general(geometry.type());
    const auto quadPos0          = referenceElement.position(0, 0);
    const auto detJ0             = geometry.integrationElement(quadPos0);

    return (transformationMatrix(geometry, quadPos0) * detJ0).eval();
  }
} // namespace Impl

/**
 * \brief Interface for displacement-based EAS elements, where linear or Green-Lagrange strains are enhanced.
 */
template <typename GEO, int ess>
struct EX
{
  static constexpr int myDim              = GEO::mydimension;
  static constexpr int strainSize         = myDim * (myDim + 1) / 2;
  static constexpr int enhancedStrainSize = ess;
  using MType                             = Eigen::Matrix<double, strainSize, enhancedStrainSize>;
  using DType                             = Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize>;

  EX() = default;
  explicit EX(const GEO& geometry)
      : geometry_{std::make_optional<GEO>(geometry)},
        T0InverseTransformed_{Impl::transformationMatrixAtCenterWithDetJ(geometry).inverse()} {}

  virtual MType calcM(const Dune::FieldVector<double, myDim>& quadPos) const = 0;

protected:
  std::optional<GEO> geometry_;
  Eigen::Matrix<double, strainSize, strainSize> T0InverseTransformed_;
};

/**
 * \brief Dummy struct for displacement-based EAS elements, i.e. 0 enhanced modes
 */
template <typename GEO>
struct E0 : EX<GEO, 0>
{
  using Base                              = EX<GEO, 0>;
  static constexpr int myDim              = Base::myDim;
  static constexpr int strainSize         = Base::strainSize;
  static constexpr int enhancedStrainSize = Base::enhancedStrainSize;
  using MType                             = typename Base::MType;
  using DType                             = typename Base::DType;

  E0() = default;
  explicit E0(const GEO& geo)
      : Base(geo) {}

  // returns an Eigen zero expression for optimization purposes
  MType calcM(const Dune::FieldVector<double, myDim>& /*quadPos*/) const override { return MType::Zero(); }
};

/**
 * \brief E4 structure for EAS with linear strains and 4 enhanced modes.
 *
 * \details The E4 structure represents an implementation of EAS for a
 * specific case where linear strains are considered, and the method includes
 * enhancement with 4 additional modes to improve the accuracy of finite element solutions.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct E4 : EX<GEO, 4>
{
  using Base                              = EX<GEO, 4>;
  static constexpr int myDim              = Base::myDim;
  static constexpr int strainSize         = Base::strainSize;
  static constexpr int enhancedStrainSize = Base::enhancedStrainSize;
  using MType                             = typename Base::MType;
  using DType                             = typename Base::DType;

  E4() = default;
  explicit E4(const GEO& geo)
      : Base(geo) {}

  MType calcM(const Dune::FieldVector<double, myDim>& quadPos) const override {
    MType M;
    M.setZero(strainSize, enhancedStrainSize);
    const double xi   = quadPos[0];
    const double eta  = quadPos[1];
    M(0, 0)           = 2 * xi - 1.0;
    M(1, 1)           = 2 * eta - 1.0;
    M(2, 2)           = 2 * xi - 1.0;
    M(2, 3)           = 2 * eta - 1.0;
    const double detJ = this->geometry_->integrationElement(quadPos);
    M                 = this->T0InverseTransformed_ / detJ * M;
    return M;
  }
};

/**
 * \brief Structure representing EAS for Q1 with 5 enhanced strains.
 *
 * This structure defines the EAS for Q1 elements with 5 enhanced strains.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct E5 : EX<GEO, 5>
{
  using Base                              = EX<GEO, 5>;
  static constexpr int myDim              = Base::myDim;
  static constexpr int strainSize         = Base::strainSize;
  static constexpr int enhancedStrainSize = Base::enhancedStrainSize;
  using MType                             = typename Base::MType;
  using DType                             = typename Base::DType;

  E5() = default;
  explicit E5(const GEO& geo)
      : Base(geo) {}

  MType calcM(const Dune::FieldVector<double, myDim>& quadPos) const override {
    MType M;
    M.setZero();
    const double xi   = quadPos[0];
    const double eta  = quadPos[1];
    M(0, 0)           = 2 * xi - 1.0;
    M(1, 1)           = 2 * eta - 1.0;
    M(2, 2)           = 2 * xi - 1.0;
    M(2, 3)           = 2 * eta - 1.0;
    M(2, 4)           = (2 * xi - 1.0) * (2 * eta - 1.0);
    const double detJ = this->geometry_->integrationElement(quadPos);
    M                 = this->T0InverseTransformed_ / detJ * M;
    return M;
  }
};

/**
 * \brief Structure representing EAS for Q1 with 7 enhanced strains.
 *
 * This structure defines the EAS for Q1 elements with 7 enhanced strains.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct E7 : EX<GEO, 7>
{
  using Base                              = EX<GEO, 7>;
  static constexpr int myDim              = Base::myDim;
  static constexpr int strainSize         = Base::strainSize;
  static constexpr int enhancedStrainSize = Base::enhancedStrainSize;
  using MType                             = typename Base::MType;
  using DType                             = typename Base::DType;

  E7() = default;
  explicit E7(const GEO& geo)
      : Base(geo) {}

  MType calcM(const Dune::FieldVector<double, myDim>& quadPos) const override {
    MType M;
    M.setZero();
    const double xi   = quadPos[0];
    const double eta  = quadPos[1];
    M(0, 0)           = 2 * xi - 1.0;
    M(1, 1)           = 2 * eta - 1.0;
    M(2, 2)           = 2 * xi - 1.0;
    M(2, 3)           = 2 * eta - 1.0;
    M(0, 4)           = (2 * xi - 1.0) * (2 * eta - 1.0);
    M(1, 5)           = (2 * xi - 1.0) * (2 * eta - 1.0);
    M(2, 6)           = (2 * xi - 1.0) * (2 * eta - 1.0);
    const double detJ = this->geometry_->integrationElement(quadPos);
    M                 = this->T0InverseTransformed_ / detJ * M;
    return M;
  }
};

/**
 * \brief Structure representing EAS for H1 with 9 enhanced strains.
 *
 * This structure defines the EAS for H1 elements with 9 enhanced strains.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct E9 : EX<GEO, 9>
{
  using Base                              = EX<GEO, 9>;
  static constexpr int myDim              = Base::myDim;
  static constexpr int strainSize         = Base::strainSize;
  static constexpr int enhancedStrainSize = Base::enhancedStrainSize;
  using MType                             = typename Base::MType;
  using DType                             = typename Base::DType;

  E9() = default;
  explicit E9(const GEO& geo)
      : Base(geo) {}

  MType calcM(const Dune::FieldVector<double, myDim>& quadPos) const override {
    MType M;
    M.setZero();
    const double xi   = quadPos[0];
    const double eta  = quadPos[1];
    const double zeta = quadPos[2];
    M(0, 0)           = 2 * xi - 1.0;
    M(1, 1)           = 2 * eta - 1.0;
    M(2, 2)           = 2 * zeta - 1.0;
    M(3, 3)           = 2 * xi - 1.0;
    M(3, 4)           = 2 * eta - 1.0;
    M(4, 5)           = 2 * xi - 1.0;
    M(4, 6)           = 2 * zeta - 1.0;
    M(5, 7)           = 2 * eta - 1.0;
    M(5, 8)           = 2 * zeta - 1.0;
    const double detJ = this->geometry_->integrationElement(quadPos);
    M                 = this->T0InverseTransformed_ / detJ * M;
    return M;
  }
};

/**
 * \brief Structure representing EAS for H1 with 21 enhanced strains.
 *
 * This structure defines the EAS for H1 elements with 21 enhanced strains.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct E21 : EX<GEO, 21>
{
  using Base                              = EX<GEO, 21>;
  static constexpr int myDim              = Base::myDim;
  static constexpr int strainSize         = Base::strainSize;
  static constexpr int enhancedStrainSize = Base::enhancedStrainSize;
  using MType                             = typename Base::MType;
  using DType                             = typename Base::DType;

  E21() = default;
  explicit E21(const GEO& geo)
      : Base(geo) {}

  MType calcM(const Dune::FieldVector<double, myDim>& quadPos) const override {
    MType M;
    M.setZero();
    const double xi   = quadPos[0];
    const double eta  = quadPos[1];
    const double zeta = quadPos[2];
    M(0, 0)           = 2 * xi - 1.0;
    M(1, 1)           = 2 * eta - 1.0;
    M(2, 2)           = 2 * zeta - 1.0;
    M(3, 3)           = 2 * xi - 1.0;
    M(3, 4)           = 2 * eta - 1.0;
    M(4, 5)           = 2 * xi - 1.0;
    M(4, 6)           = 2 * zeta - 1.0;
    M(5, 7)           = 2 * eta - 1.0;
    M(5, 8)           = 2 * zeta - 1.0;

    M(3, 9)  = (2 * xi - 1.0) * (2 * zeta - 1.0);
    M(3, 10) = (2 * eta - 1.0) * (2 * zeta - 1.0);
    M(4, 11) = (2 * xi - 1.0) * (2 * eta - 1.0);
    M(4, 12) = (2 * eta - 1.0) * (2 * zeta - 1.0);
    M(5, 13) = (2 * xi - 1.0) * (2 * eta - 1.0);
    M(5, 14) = (2 * xi - 1.0) * (2 * zeta - 1.0);

    M(0, 15) = (2 * xi - 1.0) * (2 * eta - 1.0);
    M(0, 16) = (2 * xi - 1.0) * (2 * zeta - 1.0);
    M(1, 17) = (2 * xi - 1.0) * (2 * eta - 1.0);
    M(1, 18) = (2 * eta - 1.0) * (2 * zeta - 1.0);
    M(2, 19) = (2 * xi - 1.0) * (2 * zeta - 1.0);
    M(2, 20) = (2 * eta - 1.0) * (2 * zeta - 1.0);

    const double detJ = this->geometry_->integrationElement(quadPos);
    M                 = this->T0InverseTransformed_ / detJ * M;
    return M;
  }
};

} // namespace Ikarus::EAS
