// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file polyfit.hh
 * \brief Polynomial fitting of data
 */

#pragma once

#include <dune/functions/analyticfunctions/polynomial.hh>

#include <Eigen/Dense>
namespace Ikarus::utils {
/**
 * \brief Fits a polynomial of a given degree to the given data points.
 * \ingroup utils
 * \param x The input vector of x-coordinates.
 * \param y The input vector of y-coordinates.
 * \param deg The degree of the polynomial to fit.
 * \return std::tuple<Dune::Functions::Polynomial<double>, double> A tuple containing the fitted polynomial and the
 * least square error.
 */
std::tuple<Dune::Functions::Polynomial<double>, double> polyfit(const Eigen::Ref<const Eigen::VectorXd>& x,
                                                                const Eigen::Ref<const Eigen::VectorXd>& y, int deg);
} // namespace Ikarus::utils
