// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file io.hh
 * \brief Additional python bindings for submodule io
 */
#include "../pythonhelpers.hh"

#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/pybind11.h>

#include <ikarus/io/vtkdatatag.hh>

void addBindingsToIO() {
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Ikarus;

  auto io = pybind11::module::import("ikarus.io");

  using namespace Ikarus::Vtk;
  using namespace Ikarus;
  ENUM_BINDINGS_WITH_MODULE(DataTag, io);
}