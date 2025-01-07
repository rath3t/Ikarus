// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/istl/bvector.hh>

namespace Ikarus {
template <typename Manifold>
class ValueFactory
{
  using TargetSpace = Manifold;

public:
  static void construct(Dune::BlockVector<TargetSpace>& values, const int testPointsSize = 10) {
    values.resize(testPointsSize);
    std::ranges::generate(values, []() { return TargetSpace(TargetSpace::CoordinateType::Random()); });
  }
};
} // namespace Ikarus
