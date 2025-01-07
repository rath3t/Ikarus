// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file basis.hh
 * \brief Python bindings for Ikarus basis
 */

#pragma once

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/basis.hh>

namespace Ikarus::Python {

/**
 * \brief Register a Python wrapper for an Ikarus basis class.
 *
 * \tparam BasisHandler The Ikarus handler basis class to be registered.
 * \tparam options Variadic template parameters for additional options when defining the Python class.
 *
 * \param scope A Pybind11 handle representing the Python scope where the class should be registered.
 * \param cls The Pybind11 class template to be used for registering the Ikarus basis handler class.
 *
 * \ingroup pythonbindings
 */
template <class BasisHandler, class... options>
void registerBasisHandler(pybind11::handle scope, pybind11::class_<BasisHandler, options...> cls) {
  using pybind11::operator""_a;

  using GridView       = typename BasisHandler::GridView;
  using PreBasis       = typename BasisHandler::PreBasis;
  using UntouchedBasis = typename BasisHandler::UntouchedBasis;
  using FlatBasis      = typename BasisHandler::FlatBasis;

  cls.def(pybind11::init([](const UntouchedBasis& gb) { return new BasisHandler(gb.preBasis()); }));

  cls.def("flat", [](BasisHandler& self) { return self.flat(); }, pybind11::return_value_policy::reference);

  cls.def("untouched", [](BasisHandler& self) { return self.untouched(); }, pybind11::return_value_policy::reference);
}

} // namespace Ikarus::Python
