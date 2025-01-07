// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/pybind11.h>

#include <ikarus/utils/tensorutils.hh>

void addBindingsToUtils() {
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Ikarus;

  auto utils = pybind11::module::import("ikarus.utils");

  /**
   * \brief Converts a square 1x1, 2x2 or 3x3 matrix to a Voigt notation vector.
   * \ingroup utils
   * \param E Input matrix of size (size x size).
   * \param isStrain Flag indicating whether the conversion is for strain (true) or not (false) (default is true)..
   * \return  Vector with components in Voigt notation vector.
   *
   * This function converts a square matrix to a Voigt notation vector, which contains the unique components of
   * the input matrix.
   *
   * The optional isStrain parameter allows the user to specify whether the conversion is intended for strain
   * calculations. If isStrain is true, the off-diagonal components are multiplied by 2, providing the correct Voigt
   * notation for symmetric strain tensors.
   */
  utils.def(
      "toVoigt",
      [](Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 3, 3> mat, bool isStrain = true) {
        return toVoigt(mat, isStrain);
      },
      py::arg("matrix"), py::arg("isStrain") = true);

  /**
   * \brief Converts a vector given in Voigt notation to a matrix.
   * \ingroup utils
   * \param EVoigt Voigt notation vector.
   * \param isStrain Flag indicating whether the vector represents a strain (default is true).
   * \return Matrix corresponding to the vector in Voigt notation.
   * \details
   * This function converts a vector given in Voigt notation to the corresponding matrix. The conversion depends on the
   * size The parameter `isStrain` is used to determine the conversion factor for off-diagonal components, which need to
   * be divided by 2 in the matrix representation if the quantity is a strain tensor.
   *
   * The function requires that the size of the Voigt notation vector is valid (1, 3, or 6).
   */
  utils.def(
      "fromVoigt",
      [](Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 6, 1> vec, bool isStrain = true) {
        return fromVoigt(vec, isStrain);
      },
      py::arg("vector"), py::arg("isStrain") = true);
}