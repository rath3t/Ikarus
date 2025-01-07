// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "findlinesegment.hh"

namespace Ikarus::utils {
/*
 * This function is inspired from
 * https://github.com/NicolasBoumal/manopt/blob/master/manopt/tools/identify_linear_piece.m
 */
std::tuple<Dune::Functions::Polynomial<double>, decltype(Eigen::seq(0, 0))> findLineSegment(const Eigen::VectorXd& x,
                                                                                            const Eigen::VectorXd& y,
                                                                                            int segmentSize) {
  Eigen::VectorXd errors = Eigen::VectorXd::Zero(x.size() - segmentSize);
  std::vector<Dune::Functions::Polynomial<double>> lines;
  for (int i = 0; i < errors.size(); ++i) {
    auto range = Eigen::seq(i, i + segmentSize);

    auto [poly, error] = polyfit(x(range), y(range), 1);
    errors(i)          = error;
    lines.push_back(poly);
  }
  std::forward_iterator auto minEle = std::ranges::min_element(errors.begin(), errors.end());
  auto index                        = std::distance(errors.begin(), minEle);
  auto range                        = Eigen::seq(index, index + segmentSize);
  return std::make_tuple(lines[index], range);
}
} // namespace Ikarus::utils
