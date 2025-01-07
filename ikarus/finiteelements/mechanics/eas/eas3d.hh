// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file eas3d.hh
 * \brief Definition of the types of EAS formulations for 3D elements.
 * \ingroup  mechanics
 */

#pragma once

#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::EAS {

/**
 * \brief Structure representing EAS for H1 with 9 enhanced strains.
 *
 * This structure defines the EAS for H1 elements with 9 enhanced strains.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct H1E9
{
  static constexpr int strainSize         = 6;
  static constexpr int enhancedStrainSize = 9;
  using MType                             = Eigen::Matrix<double, strainSize, enhancedStrainSize>;
  using DType                             = Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize>;

  H1E9() = default;
  explicit H1E9(const GEO& geometry)
      : geometry_{std::make_optional<GEO>(geometry)},
        T0InverseTransformed_{calcTransformationMatrix3D(geometry)} {}

  auto calcM(const Dune::FieldVector<double, 3>& quadPos) const {
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
    const double detJ = geometry_->integrationElement(quadPos);
    M                 = T0InverseTransformed_ / detJ * M;
    return M;
  }

private:
  std::optional<GEO> geometry_;
  Eigen::Matrix<double, strainSize, strainSize> T0InverseTransformed_;
};

/**
 * \brief Structure representing EAS for H1 with 21 enhanced strains.
 *
 * This structure defines the EAS for H1 elements with 21 enhanced strains.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct H1E21
{
  static constexpr int strainSize         = 6;
  static constexpr int enhancedStrainSize = 21;
  using MType                             = Eigen::Matrix<double, strainSize, enhancedStrainSize>;
  using DType                             = Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize>;

  H1E21() = default;
  explicit H1E21(const GEO& geometry)
      : geometry_{std::make_optional<GEO>(geometry)},
        T0InverseTransformed{calcTransformationMatrix3D(geometry)} {}

  auto calcM(const Dune::FieldVector<double, 3>& quadPos) const {
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

    const double detJ = geometry_->integrationElement(quadPos);
    M                 = T0InverseTransformed / detJ * M;
    return M;
  }

private:
  std::optional<GEO> geometry_;
  Eigen::Matrix<double, strainSize, strainSize> T0InverseTransformed;
};
} // namespace Ikarus::EAS
