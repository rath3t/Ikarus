// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <muesli/muesli.h>

#include <dune/common/hybridutilities.hh>

#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Concepts {
template <typename MAT>
concept MuesliMaterialImplementation = requires { std::is_base_of_v<muesli::material, MAT>; };
} // namespace Ikarus::Concepts

namespace Ikarus::Materials::Muesli {

// Alias for Muesli material properties
using MaterialProperties = muesli::materialProperties;

template <MPTuple MPT>
inline MaterialProperties propertiesFromIkarusMaterialParameters(const MPT& mpt) {
  auto converter = convertLameConstants(mpt);

  auto mpm = muesli::materialProperties{};
  mpm.insert({"lambda", converter.toLamesFirstParameter()});
  mpm.insert({"mu", converter.toShearModulus()});

  return mpm;
}

inline void addRegularizedTag(MaterialProperties& mpm) { mpm.insert({"subtype regularized", 0}); }
inline void addCompressibleTag(MaterialProperties& mpm) { mpm.insert({"compressible", 0}); }
inline void addIncompressibleTag(MaterialProperties& mpm) { mpm.insert({"incompressible", 0}); }

template <typename Derived>
inline void toistensor(istensor& it, const Eigen::MatrixBase<Derived>& C) {
  it = istensor(C(0, 0), C(1, 1), C(2, 2), C(1, 2), C(2, 0), C(0, 1));
}

template <typename ScalarType = double, int dim = 3>
inline Eigen::Matrix<ScalarType, dim, dim> toMatrix(const istensor& it) {
  auto S = Eigen::Matrix<double, dim, dim>{};
  for (auto i : Dune::Hybrid::integralRange(dim))
    for (auto j : Dune::Hybrid::integralRange(dim))
      S(i, j) = it(i, j);
  return S;
}

template <typename ScalarType = double, std::size_t dim = 3>
inline auto toTensor(const itensor4& it) -> Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>> {
  Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>> moduli{};
  moduli.setZero();
  for (auto i : Dune::range(3))
    for (auto j : Dune::range(3))
      for (auto k : Dune::range(3))
        for (auto l : Dune::range(3))
          moduli(i, j, k, l) = it(i, j, k, l);

  return moduli;
}

template <Concepts::MuesliMaterialImplementation MAT>
constexpr std::string materialName() {
  std::string matName = Dune::className<MAT>();

  // Find and remove the "muesli::" namespace
  std::string prefix = "muesli::";
  matName            = matName.substr(prefix.size());

  // Capitalize the first letter
  if (!matName.empty() && std::islower(matName[0])) {
    matName[0] = std::toupper(matName[0]);
  }
  return matName;
}

} // namespace Ikarus::Materials::Muesli