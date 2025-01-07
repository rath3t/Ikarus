// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file pythonautodiffdefinitions.hh
 * \brief Implementation of forwarding autodiff types from python to c++ and vice versa
 */

#pragma once
#include <cstddef>
#include <Python.h>

#include <dune/fufem/dunepython.hh>

#include <autodiff/forward/real/real.hpp>

namespace Python {
// *****************************************************************************
// specializations of Conversion that do the PyObject* <-> C++-type conversion
// *****************************************************************************

/**
 * \brief Conversion specialization for autodiff::Real type.
 * \tparam order The order of autodiff::Real.
 * \tparam T The underlying type of autodiff::Real.
 */
template <std::size_t order, class T>
struct Conversion<autodiff::Real<order, T>>
{
  enum
  {
    useDefaultConstructorConversion = true
  };

  /**
   * \brief Convert autodiff::Real to PyObject*.
   * \param list Python list object to be populated.
   * \param v The autodiff::Real object to be converted.
   */
  static void toC(PyObject* list, autodiff::Real<order, T>& v) {
    auto rlist = Reference(Imp::inc(list));
    // Reference is needed to enable the ".get()" function and "Imp::inc" is
    // needed since Reference owns the PyObject pointer and decrements it at the end of the scope
    // Imp::inc artificially increases the reference counter by one.
    // When we return from this function, these two cancel out and the PyObject* is as before

    auto wF = Callable(rlist.get("__getitem__"));
    for (std::size_t i = 0; i < order + 1; ++i)
      v[i] = PyFloat_AsDouble(wF(i));
  }

  /**
   * \brief Convert PyObject* to autodiff::Real.
   * \param v The autodiff::Real object to be populated.
   * \return PyObject* representing the autodiff::Real object.
   */
  static PyObject* toPy(const autodiff::Real<order, T>& v) {
    auto pyMain           = Python::main();
    Python::Module module = pyMain.import("autodiff");

    auto real1stClass = module.get("real1st");
    auto real1st      = Callable(Imp::inc(real1stClass))();
    auto wF           = Callable(Imp::inc(real1st).get("__setitem__"));
    for (std::size_t i = 0; i < order + 1; ++i)
      wF(i, v[i]); // real1st.__setitem__(i,v[i]);

    return real1st;
  }
};
} // namespace Python
