// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file eas2d.hh
 * \brief Definition of the types of EAS formulations for 2D elements.
 * \ingroup  eas
 */

#pragma once

#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::EAS {

/**
 * \brief Dummy struct for displacement-based EAS elements, i.e. 0 enhanced modes
 */
template <typename GEO>
struct E0
{
  static constexpr int strainSize         = GEO::mydimension * (GEO::mydimension + 1) / 2;
  static constexpr int enhancedStrainSize = 0;
  using MType                             = Eigen::Matrix<double, strainSize, enhancedStrainSize>;
  using DType                             = Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize>;

  E0() = default;
  explicit E0(const GEO& /*geometry*/) {}

  // returns an Eigen zero expression for optimization purposes
  static auto calcM(const Dune::FieldVector<double, GEO::mydimension>& /*quadPos*/) { return MType::Zero(); }
};

/**
 * \brief Q1E4 structure for EAS with linear strains and 4 enhanced modes.
 *
 * \details The Q1E4 structure represents an implementation of EAS for a
 * specific case where linear strains are considered, and the method includes
 * enhancement with 4 additional modes to improve the accuracy of finite element solutions.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct Q1E4
{
  static constexpr int strainSize         = 3;
  static constexpr int enhancedStrainSize = 4;
  using MType                             = Eigen::Matrix<double, strainSize, enhancedStrainSize>;
  using DType                             = Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize>;

  Q1E4() = default;
  explicit Q1E4(const GEO& geometry)
      : geometry_{std::make_optional<GEO>(geometry)},
        T0InverseTransformed_{calcTransformationMatrix2D(geometry)} {}

  auto calcM(const Dune::FieldVector<double, 2>& quadPos) const {
    MType M;
    M.setZero(strainSize, enhancedStrainSize);
    const double xi   = quadPos[0];
    const double eta  = quadPos[1];
    M(0, 0)           = 2 * xi - 1.0;
    M(1, 1)           = 2 * eta - 1.0;
    M(2, 2)           = 2 * xi - 1.0;
    M(2, 3)           = 2 * eta - 1.0;
    const double detJ = geometry_->integrationElement(quadPos);
    M                 = T0InverseTransformed_ / detJ * M;
    return M;
  }

private:
  std::optional<GEO> geometry_;
  Eigen::Matrix3d T0InverseTransformed_;
};

/**
 * \brief Structure representing EAS for Q1 with 5 enhanced strains.
 *
 * This structure defines the EAS for Q1 elements with 5 enhanced strains.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct Q1E5
{
  static constexpr int strainSize         = 3;
  static constexpr int enhancedStrainSize = 5;
  using MType                             = Eigen::Matrix<double, strainSize, enhancedStrainSize>;
  using DType                             = Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize>;

  Q1E5() = default;
  explicit Q1E5(const GEO& geometry)
      : geometry_{std::make_optional<GEO>(geometry)},
        T0InverseTransformed_{calcTransformationMatrix2D(geometry)} {}

  auto calcM(const Dune::FieldVector<double, 2>& quadPos) const {
    MType M;
    M.setZero();
    const double xi   = quadPos[0];
    const double eta  = quadPos[1];
    M(0, 0)           = 2 * xi - 1.0;
    M(1, 1)           = 2 * eta - 1.0;
    M(2, 2)           = 2 * xi - 1.0;
    M(2, 3)           = 2 * eta - 1.0;
    M(2, 4)           = (2 * xi - 1.0) * (2 * eta - 1.0);
    const double detJ = geometry_->integrationElement(quadPos);
    M                 = T0InverseTransformed_ / detJ * M;
    return M;
  }

private:
  std::optional<GEO> geometry_;
  Eigen::Matrix3d T0InverseTransformed_;
};

/**
 * \brief Structure representing EAS for Q1 with 7 enhanced strains.
 *
 * This structure defines the EAS for Q1 elements with 7 enhanced strains.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct Q1E7
{
  static constexpr int strainSize         = 3;
  static constexpr int enhancedStrainSize = 7;
  using MType                             = Eigen::Matrix<double, strainSize, enhancedStrainSize>;
  using DType                             = Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize>;

  Q1E7() = default;
  explicit Q1E7(const GEO& geometry)
      : geometry_{std::make_optional<GEO>(geometry)},
        T0InverseTransformed_{calcTransformationMatrix2D(geometry)} {}

  auto calcM(const Dune::FieldVector<double, 2>& quadPos) const {
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
    const double detJ = geometry_->integrationElement(quadPos);
    M                 = T0InverseTransformed_ / detJ * M;
    return M;
  }

private:
  std::optional<GEO> geometry_;
  Eigen::Matrix3d T0InverseTransformed_;
};

} // namespace Ikarus::EAS
