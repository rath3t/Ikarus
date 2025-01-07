// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/common/bitsetvector.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <ikarus/utils/modalanalysis/lumpingschemes.hh>

namespace Ikarus::Python {

template <class ModalAnalysis, class... options>
void registerModalAnalysis(pybind11::handle scope, pybind11::class_<ModalAnalysis, options...> cls) {
  using pybind11::operator""_a;

  using Assembler           = typename ModalAnalysis::Assembler;
  using DirichletValuesType = typename Assembler::DirichletValuesType;
  using FEContainer         = typename ModalAnalysis::FEContainer;

  cls.def(pybind11::init([](const pybind11::list& fes, const DirichletValuesType& dirichletValues) {
            FEContainer fesV = fes.template cast<FEContainer>();
            return new ModalAnalysis(std::move(fesV), dirichletValues);
          }),
          pybind11::keep_alive<1, 3>());

  cls.def("compute", &ModalAnalysis::compute, pybind11::return_value_policy::copy, pybind11::arg("maxit") = 1000,
          pybind11::arg("tolerance") = 1e-10);
  cls.def("angularFrequencies", &ModalAnalysis::angularFrequencies, pybind11::return_value_policy::copy);
  cls.def("naturalFrequencies", &ModalAnalysis::naturalFrequencies, pybind11::return_value_policy::copy);
  cls.def("squaredAngularFrequencies", &ModalAnalysis::squaredAngularFrequencies, pybind11::return_value_policy::copy);
  cls.def("eigenmodes", &ModalAnalysis::eigenmodes, pybind11::return_value_policy::copy);
  cls.def_property_readonly("nev", &ModalAnalysis::nev);
  cls.def("writeEigenModes", &ModalAnalysis::writeEigenModes, pybind11::arg("filename"),
          pybind11::arg("nev") = std::nullopt);

  // For now, only row-sum-lumping is possible through the python bindings.
  cls.def("bindLumpingScheme",
          [](ModalAnalysis& self) { self.template bindLumpingScheme<Dynamics::LumpingSchemes::RowSumLumping>(); });
  cls.def("unBindLumpingScheme", &ModalAnalysis::unBindLumpingSchemes);
}

} // namespace Ikarus::Python
