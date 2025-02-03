// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "io/io.hh"
#include "materials/materials.hh"
#include "pythonhelpers.hh"
#include "utils/utils.hh"

#include <dune/common/float_cmp.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/pybind11.h>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/python/finiteelements/scalarwrapper.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>

/**
 * \brief Registers the ScalarWrapper class template with pybind11, adding various operations and constructors.
 *
 * This function template registers the ScalarWrapper class in the given Python module with specified names for
 * the class and its scalar type. It adds support for basic arithmetic operations, comparisons, and string
 * representation.
 *
 * \tparam ScalarType The type of the scalar value to be wrapped.
 * \param m The pybind11 module where the class will be registered.
 * \param name The name of the class in the Python module.
 * \param typeName The cpp name of the scalar type in the Python module.
 */
template <typename ScalarType>
void registerScalarWrapper(pybind11::module& m, std::string name, std::string typeName) {
  namespace py = pybind11;
  using namespace Ikarus;
  using namespace Eigen;
  using VWd           = ScalarWrapper<ScalarType>;
  using RawScalarType = Dune::ResolveRef_t<ScalarType>;

  auto scalarWrapper = Dune::Python::insertClass<VWd>(m, name, Dune::Python::GenerateTypeName(typeName)).first;
  if constexpr (std::is_same_v<VWd, RawScalarType>)
    scalarWrapper.def(py::init<RawScalarType>(), "Constructor with a raw scalar value.");
  scalarWrapper.def(py::init<VWd>(), "Copy constructor.");
  scalarWrapper.def(py::init<ScalarType>(), "Constructor with a scalar value.");
  scalarWrapper.def("__repr__", [](const VWd& d) { return std::to_string(d.value()); }, "String representation.");
  scalarWrapper.def(
      "__eq__", [](const VWd& x, const VWd& y) { return Dune::FloatCmp::eq(x.value(), y.value()); },
      "Equality comparison with another ScalarWrapper.");
  scalarWrapper.def(
      "__eq__", [](const VWd& x, const RawScalarType& y) { return Dune::FloatCmp::eq(x.value(), y); },
      "Equality comparison with a raw scalar value.");
  scalarWrapper.def(
      "__eq__", [](const RawScalarType& x, const VWd& y) { return Dune::FloatCmp::eq(x, y.value()); },
      "Equality comparison of a raw scalar value with a ScalarWrapper.");
  scalarWrapper.def(py::self + py::self, "Addition of two ScalarWrapper instances.")
      .def(py::self - py::self, "Subtraction of two ScalarWrapper instances.")
      .def(py::self += py::self, "Addition assignment of another ScalarWrapper instance.")
      .def(py::self -= py::self, "Subtraction assignment of another ScalarWrapper instance.")
      .def(-py::self, "Negation of the ScalarWrapper value.")
      .def(py::self *= RawScalarType(), "Multiplication assignment with a raw scalar value.")
      .def(py::self /= RawScalarType(), "Division assignment with a raw scalar value.")
      .def(RawScalarType() * py::self, "Multiplication of a raw scalar value with a ScalarWrapper.")
      .def(py::self * RawScalarType(), "Multiplication of a ScalarWrapper with a raw scalar value.")
      .doc() =
      R"(Since Python does not support passing python float by reference to a double&, we have to wrap everything to allow mutable scalar types that can be passed back and forth.
              See also: https://pybind11.readthedocs.io/en/stable/faq.html#limitations-involving-reference-arguments)";
}

PYBIND11_MODULE(_ikarus, m) {
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Ikarus;
  using namespace Eigen;

  ENUM_BINDINGS(ScalarAffordance);
  ENUM_BINDINGS(VectorAffordance);
  ENUM_BINDINGS(MatrixAffordance);
  ENUM_BINDINGS(FESolutions);
  ENUM_BINDINGS(FEParameter);
  ENUM_BINDINGS(SolverTypeTag);
  ENUM_BINDINGS(DBCOption);

  py::class_<std::tuple<ScalarAffordance, VectorAffordance, MatrixAffordance>>(m, "Base");
  auto affordanceCollections3 =
      Dune::Python::insertClass<AffordanceCollection<ScalarAffordance, VectorAffordance, MatrixAffordance>,
                                std::tuple<ScalarAffordance, VectorAffordance, MatrixAffordance>>(
          m, "AffordanceCollection",
          Dune::Python::GenerateTypeName("AffordanceCollection<ScalarAffordance,VectorAffordance,MatrixAffordance>"),
          Dune::Python::IncludeFiles("ikarus/finiteelements/ferequirements.hh"))
          .first;

  affordanceCollections3.def(py::init<ScalarAffordance, VectorAffordance, MatrixAffordance>());
  affordanceCollections3.def_property_readonly_static("elastoStatics",
                                                      [](py::object) { return AffordanceCollections::elastoStatics; });

  registerScalarWrapper<double>(m, "Scalar", "ScalarWrapper<double>");
  registerScalarWrapper<std::reference_wrapper<double>>(m, "ScalarRef",
                                                        "ScalarWrapper<std::reference_wrapper<double>>");

  addBindingsToMaterials();

  addBindingsToUtils();
  addBindingsToIO();
}
