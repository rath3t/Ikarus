// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#if ENABLE_MUESLI

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

/**
 * \brief converts Ikarus material parameters to Muesli material properties
 *
 * \tparam MPT the type of the Ikarus material parameters
 * \param mpt the Ikarus material parameters
 * \return MaterialProperties meusli material properties with the Lambd and mu set
 */
template <Concepts::MPTuple MPT>
inline MaterialProperties propertiesFromIkarusMaterialParameters(const MPT& mpt) {
  auto converter = convertLameConstants(mpt);

  auto mpm = muesli::materialProperties{};
  mpm.insert({"lambda", converter.toLamesFirstParameter()});
  mpm.insert({"mu", converter.toShearModulus()});

  return mpm;
}

/**
 * \brief adds a specific tag to the muesli materialproperties, thic can be used for example to add the `regularized`
 * tag for Neo-Hooke or the `compressible` tag for Yeoh and Arruda-Boyce
 */
inline void addTag(MaterialProperties& mpm, const std::string& tagName, double tagValue = 0) {
  mpm.insert({tagName, tagValue});
}

/**
 * \brief Converts the entries of a Eigen::Matrix to a provided muesli::tensor (symmetric 2nd order tensor)
 *
 * \tparam Derived the derived Eigen::Matrix type
 * \param it provided istensor
 * \param C the Eigen::Matrix that is to be converted
 */
template <typename Derived>
inline void toistensor(istensor& it, const Eigen::MatrixBase<Derived>& C) {
  it = istensor(C(0, 0), C(1, 1), C(2, 2), C(1, 2), C(2, 0), C(0, 1));
}

/**
 * \brief Converts a provided muesli::istensor (symmetric 2nd order tensor) to a Eigen::Matrix
 *
 * \tparam ScalarType the ScalarType (defaults to double)
 * \tparam dim the dimension (defaults to 3)
 * \param it provided istensor
 * \return Eigen::Matrix<ScalarType, dim, dim> converted matrix
 */
template <typename ScalarType = double, int dim = 3>
inline Eigen::Matrix<ScalarType, dim, dim> toMatrix(const istensor& it) {
  auto S = Eigen::Matrix<double, dim, dim>{};
  for (auto i : Dune::Hybrid::integralRange(dim))
    for (auto j : Dune::Hybrid::integralRange(dim))
      S(i, j) = it(i, j);
  return S;
}

/**
 * \brief Converts a provided muesli::istensor4 (symmetric 4th order tensor) to a Eigen::TensorFixedSize
 *
 * \tparam ScalarType the ScalaType (defaults to double)
 * \tparam dim the dimension (defaults to 3)
 * \param it provided istensor
 * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>> converted tensor
 */
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

#else
  #error Muesli materials depends on the Muesli library, which is not included
#endif