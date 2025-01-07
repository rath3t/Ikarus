// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/python/pybind11/pybind11.h>

namespace Ikarus::Python {

/**
 * @brief Registers a NonLinearElasticPre class in Python.
 *
 * @tparam NonLinearElasticPre The NonLinearElasticPre class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class NonLinearElasticPre, class... options>
void registerNonLinearElasticPre(pybind11::handle scope, pybind11::class_<NonLinearElasticPre, options...> cls) {
  using Material = typename NonLinearElasticPre::Material;
  cls.def(pybind11::init([](const Material& mat, double density) { return new NonLinearElasticPre(mat, density); }));
}

/**
 * @brief Registers a LinearElasticPre class in Python.
 *
 * @tparam LinearElasticPre The LinearElasticPre class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class LinearElasticPre, class... options>
void registerLinearElasticPre(pybind11::handle scope, pybind11::class_<LinearElasticPre, options...> cls) {
  using Material = typename LinearElasticPre::Material;
  cls.def(pybind11::init([](const Material& mat, double density) { return new LinearElasticPre(mat, density); }));
}

/**
 * @brief Registers a TrussPre class in Python.
 *
 * @tparam TrussPre The TrussPre class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class TrussPre, class... options>
void registerTrussPre(pybind11::handle scope, pybind11::class_<TrussPre, options...> cls) {
  cls.def(pybind11::init([](double emod, double area, double density) { return new TrussPre({emod, area, density}); }));
}

/**
 * @brief Registers a KirchhoffLoveShellPre class in Python.
 *
 * @tparam KirchhoffLoveShellPre The KirchhoffLoveShellPre class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class KirchhoffLoveShellPre, class... options>
void registerKirchhoffLoveShellPre(pybind11::handle scope, pybind11::class_<KirchhoffLoveShellPre, options...> cls) {
  cls.def(pybind11::init(
      [](double E, double nu, double h, double density) { return new KirchhoffLoveShellPre({E, nu}, h, density); }));
}

/**
 * @brief Registers an EnhancedAssumedStrainsPre class in Python.
 *
 * @tparam EASPre The EnhancedAssumedStrainsPre class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class EASPre, class... options>
void registerEnhancedAssumedStrainsPre(pybind11::handle scope, pybind11::class_<EASPre, options...> cls) {
  cls.def(pybind11::init([](int numberOfParameter) { return new EASPre(numberOfParameter); }));
}

/**
 * @brief Registers a NeumannBoundaryLoadPre class in Python.
 *
 * @tparam NeumannBoundaryLoadPre The NeumannBoundaryLoadPre class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class NeumannBoundaryLoadPre, class... options>
void registerNeumannBoundaryLoadPre(pybind11::handle scope, pybind11::class_<NeumannBoundaryLoadPre, options...> cls) {
  using BoundaryPatchType = typename NeumannBoundaryLoadPre::BoundaryPatchType;
  using GridView          = typename NeumannBoundaryLoadPre::GridView;

  using LoadFunction = std::function<Eigen::Vector<double, NeumannBoundaryLoadPre::worldDim>(
      Dune::FieldVector<double, NeumannBoundaryLoadPre::worldDim>, const double&)>;
  cls.def(pybind11::init([](const BoundaryPatchType& patch, LoadFunction volumeLoad) {
            return new NeumannBoundaryLoadPre(&patch, volumeLoad);
          }),
          pybind11::keep_alive<1, 2>());
}

/**
 * @brief Registers a VolumeLoadPre class in Python.
 *
 * @tparam VolumeLoadPre The VolumeLoadPre class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class VolumeLoadPre, class... options>
void registerVolumeLoadPre(pybind11::handle scope, pybind11::class_<VolumeLoadPre, options...> cls) {
  using LoadFunction = std::function<Eigen::Vector<double, VolumeLoadPre::worldDim>(
      Dune::FieldVector<double, VolumeLoadPre::worldDim>, const double&)>;
  cls.def(pybind11::init([](LoadFunction volumeLoad) { return new VolumeLoadPre(volumeLoad); }));
}

} // namespace Ikarus::Python
