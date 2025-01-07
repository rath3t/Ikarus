// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file tensorproductquadrule.hh
 * \brief Free function for creating a tensor product quadrature rule and other Dune::Quadrature rule utils
 */
#pragma once

#include <dune/geometry/quadraturerules.hh>
namespace Ikarus {

/**
 * \brief Computes the tensor product quadrature rule using a base and one-dimensional quadrature rule.
 *
 * \details This function constructs a tensor product quadrature rule by combining a base quadrature rule
 * and a one-dimensional quadrature rule. It iterates over the points of the base quadrature and
 * combines them with each point of the one-dimensional quadrature.
 *
 * \tparam BaseQuadrature Type of the base quadrature.
 * \tparam Quadrature Type of the one-dimensional quadrature.
 *
 * \param baseQuad The base quadrature rule.
 * \param onedQuad The one-dimensional quadrature rule.
 *
 * \return Dune::QuadratureRule<double, BaseQuadrature::d + 1> - Tensor product quadrature rule.
 */
template <class BaseQuadrature, class Quadrature>
auto tensorProductQuadrature(const BaseQuadrature& baseQuad, const Quadrature& onedQuad) {
  constexpr int baseQuadDim       = BaseQuadrature::d;
  auto rule                       = Dune::QuadratureRule<double, baseQuadDim + 1>();
  const unsigned int baseQuadSize = baseQuad.size();
  for (unsigned int bqi = 0; bqi < baseQuadSize; ++bqi) {
    const typename Dune::QuadraturePoint<double, baseQuadDim>::Vector& basePoint = baseQuad[bqi].position();
    const double& baseWeight                                                     = baseQuad[bqi].weight();

    typename Dune::QuadraturePoint<double, baseQuadDim + 1>::Vector point;
    for (unsigned int i = 0; i < baseQuadDim; ++i)
      point[i] = basePoint[i];

    const unsigned int onedQuadSize = onedQuad.size();
    for (unsigned int oqi = 0; oqi < onedQuadSize; ++oqi) {
      point[baseQuadDim] = onedQuad[oqi].position()[0];
      rule.emplace_back(Dune::QuadraturePoint(point, baseWeight * onedQuad[oqi].weight()));
    }
  }
  return rule;
}
} // namespace Ikarus
