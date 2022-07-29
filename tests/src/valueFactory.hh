//
// Created by lex on 7/29/22.
//

#pragma once

#include <dune/istl/bvector.hh>

namespace Ikarus {
  template <typename Manifold>
  class ValueFactory
  {
    using TargetSpace = Manifold;

  public:
    static void construct(Dune::BlockVector<TargetSpace >& values) {
      const int  testPointsSize=10;
      std::ranges::generate_n(values.begin(),testPointsSize,[](){return TargetSpace(TargetSpace::CoordinateType::Random() );});
    }
  };
}  // namespace Ikarus
