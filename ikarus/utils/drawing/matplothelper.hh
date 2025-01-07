// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file matplothelper.hh
 * \brief Helper function for the matplot library
 * \ingroup io
 */

#pragma once
#include <matplot/matplot.h>

#include <Eigen/Core>
namespace Ikarus::plot {
/**
 * \brief Draw a 2D plot with given x and y vectors.
 *
 * This function uses the Matplot library to create a 2D plot with the provided x and y vectors.
 *
 * \param x Vector representing the x-axis values.
 * \param y Vector representing the y-axis values.
 */
void draw_xy(const Eigen::VectorXd& x, const Eigen::VectorXd& y);

/**
 * \brief Draw the plot of a mathematical function.
 *
 * This function uses the Matplot library to draw the plot of a mathematical function within the specified x-range.
 *
 * \tparam FunctionType The type of the mathematical function.
 * \param f The mathematical function to be plotted.
 * \param xRange A pair representing the range of x-axis values for plotting.
 * \param eValuationPoints The number of points to evaluate the function within the given range.
 */
template <typename FunctionType>
void drawFunction(FunctionType&& f, std::pair<double, double>&& xRange, int eValuationPoints = 100) {
  std::vector<double> x = matplot::linspace(xRange.first, xRange.second, eValuationPoints);
  std::vector<double> y = matplot::transform(x, [&f](auto xL) { return f(xL); });
  matplot::plot(x, y, "-o");
  matplot::hold(matplot::on);
}

} // namespace Ikarus::plot
