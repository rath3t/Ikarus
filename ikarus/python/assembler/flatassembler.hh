// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file flatassembler.hh
 * \brief Python bindings for assemblers
 */

#pragma once

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>
#include <dune/python/pybind11/stl_bind.h>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/finiteelements/feconfiguration.hh>
#include <ikarus/utils/basis.hh>

namespace Ikarus::Python {

/**
 \brief Register Python bindings for a assembler class.
\n

 The registered class will have an initializer that takes a list of finite elements (`fes`) and a
`DirichletValuesType` object. \n It exposes several member functions to Python:  \n <ul> <li> `matrix(req)`:
Returns a dense matrix based on the specified `FERequirementType`.    \n <li> `vector(req,affordance,dbcOption)`:
Returns a vector based on the specified `FERequirementType`.   \n <li> `scalar(req,affordance)`: Returns a scalar
based on the specified `FERequirementType`.    \n <li> `createFullVector(redVec)`: Creates a full vector from a
reduced vector. \n <li> `reducedSize()`: Returns the size of the reduced space.   \n
</ul> \n
 \tparam Assembler The assembler class to be registered.
 \tparam options Variadic template parameters for additional options when defining the Python class.
 \param scope A Pybind11 handle representing the Python scope where the class should be registered.
 \param cls The Pybind11 class template to be used for registering the assembler class.
 \ingroup pythonbindings
*/
template <class Assembler, class... options>
void registerFlatAssembler(pybind11::handle scope, pybind11::class_<Assembler, options...> cls) {
  using pybind11::operator""_a;
  using FEContainer              = typename Assembler::FEContainer;
  using Basis                    = typename Assembler::Basis;
  using DirichletValuesType      = typename Assembler::DirichletValuesType;
  using AffordanceCollectionType = typename Assembler::AffordanceCollectionType;
  using FERequirementType        = typename Assembler::FERequirement;
  using SizeType                 = typename Assembler::SizeType;

  cls.def(pybind11::init([](const pybind11::list& fes, const DirichletValuesType& dirichletValues) {
            /*here a copy of the whole vector of fes takes place! There is no way to prevent this if we want that
             * the user can pass native python lists here, see
             * https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html */
            FEContainer fesV = fes.template cast<FEContainer>();
            return new Assembler(std::move(fesV), dirichletValues);
          }),
          pybind11::keep_alive<1, 3>());

  /* sparse matrices need to be copied to python therefore we remove the reference of the return type, see */
  /*  https://github.com/pybind/pybind11/blob/cbb876cc7b02c5f57e715cbc2c46ead3d1fbcf79/tests/test_eigen_matrix.cpp#L332-L341
   */
  cls.def(
      "matrix",
      [](Assembler& self, const FERequirementType& req, Ikarus::MatrixAffordance affordance,
         Ikarus::DBCOption dbcOption) -> std::remove_cvref_t<decltype(self.matrix(req, affordance))> {
        return self.matrix(req, affordance, dbcOption);
      },
      pybind11::return_value_policy::copy);

  cls.def(
      "matrix", [](Assembler& self) -> std::remove_cvref_t<decltype(self.matrix())> { return self.matrix(); },
      pybind11::return_value_policy::copy);

  cls.def(
      "matrix",
      [](Assembler& self, Ikarus::DBCOption dbcOption) -> std::remove_cvref_t<decltype(self.matrix(dbcOption))> {
        return self.matrix(dbcOption);
      },
      pybind11::return_value_policy::copy);

  cls.def(
      "vector",
      [](Assembler& self, const FERequirementType& req, Ikarus::VectorAffordance affordance,
         Ikarus::DBCOption dbcOption) { return self.vector(req, affordance, dbcOption); },
      pybind11::return_value_policy::reference);

  cls.def("vector", [](Assembler& self) { return self.vector(); }, pybind11::return_value_policy::reference);

  cls.def(
      "vector", [](Assembler& self, Ikarus::DBCOption dbcOption) { return self.vector(dbcOption); },
      pybind11::return_value_policy::reference);

  cls.def(
      "scalar",
      [](Assembler& self, const FERequirementType& req, Ikarus::ScalarAffordance affordance) {
        return self.scalar(req, affordance);
      },
      pybind11::return_value_policy::copy);

  cls.def("scalar", [](Assembler& self) { return self.scalar(); }, pybind11::return_value_policy::copy);

  cls.def(
      "createFullVector",
      [](Assembler& self, Eigen::Ref<const Eigen::VectorXd> redVec) { return self.createFullVector(redVec); },
      pybind11::return_value_policy::move);
  cls.def("reducedSize", [](Assembler& self) { return self.reducedSize(); }, pybind11::return_value_policy::copy);
  cls.def("bind", [](Assembler& self, const FERequirementType& req, AffordanceCollectionType affordance,
                     DBCOption dbcOption = DBCOption::Full) { return self.bind(req, affordance, dbcOption); });
  cls.def("bind", [](Assembler& self, const FERequirementType& req) { return self.bind(req); });
  cls.def("bind", [](Assembler& self, const AffordanceCollectionType affordance) { return self.bind(affordance); });
  cls.def("bind", [](Assembler& self, const DBCOption dbcOption) { return self.bind(dbcOption); });
  cls.def("bound", [](Assembler& self) { return self.bound(); });
  cls.def("boundToRequirement", [](Assembler& self) { return self.boundToRequirement(); });
  cls.def("boundToAffordanceCollection", [](Assembler& self) { return self.boundToAffordanceCollection(); });
  cls.def("boundToDBCOption", [](Assembler& self) { return self.boundToDBCOption(); });
  cls.def("requirement", [](Assembler& self) { return self.requirement(); });
  cls.def("affordanceCollection", [](Assembler& self) { return self.affordanceCollection(); });
  cls.def("dBCOption", [](Assembler& self) { return self.dBCOption(); });
  cls.def_property_readonly("size", [](Assembler& self) { return self.size(); });
  cls.def("__len__", [](Assembler& self) { return self.size(); });
  cls.def("constraintsBelow", [](Assembler& self, SizeType i) { return self.constraintsBelow(i); });
  cls.def("isConstrained", [](Assembler& self, SizeType i) { return self.isConstrained(i); });
  cls.def_property_readonly("gridView", [](Assembler& self) { return self.gridView(); });
}

#define MAKE_ASSEMBLER_REGISTERY_FUNCTION(name)                                              \
  template <class Assembler, class... options>                                               \
  void register##name(pybind11::handle scope, pybind11::class_<Assembler, options...> cls) { \
    registerFlatAssembler(scope, cls);                                                       \
  }

MAKE_ASSEMBLER_REGISTERY_FUNCTION(SparseFlatAssembler);
MAKE_ASSEMBLER_REGISTERY_FUNCTION(DenseFlatAssembler);

} // namespace Ikarus::Python
