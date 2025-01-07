// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file findlinesegment.hh
 * \brief Implementation of the findLineSegment algorithm
 */

#pragma once
#include "polyfit.hh"

#include <Eigen/Core>

namespace Ikarus::utils {
/**
 * \brief Find a linear segment in a set of data points.
 *\ingroup utils
 * \details his function is inspired by the MATLAB code at:
 * https://github.com/NicolasBoumal/manopt/blob/master/manopt/tools/identify_linear_piece.m
 * It is designed to find the most linear segment in a set of data points
 *
 * \param x The x-coordinates of the data points.
 * \param y The y-coordinates of the data points.
 * \param segmentSize The size of the line segment to be identified.
 * \return A tuple containing the polynomial representing the identified line segment and the indices of the data
 *points in the segment.
 */
std::tuple<Dune::Functions::Polynomial<double>, decltype(Eigen::seq(0, 0))> findLineSegment(const Eigen::VectorXd& x,
                                                                                            const Eigen::VectorXd& y,
                                                                                            int segmentSize);
} // namespace Ikarus::utils
