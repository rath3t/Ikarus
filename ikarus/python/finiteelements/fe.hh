// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file traction.hh
 * \brief Python bindings for the traction pre
 */

#pragma once

#include <dune/python/common/typeregistry.hh>
#include <dune/python/functions/globalbasis.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/python/finiteelements/registerferequirements.hh>

namespace Ikarus::Python {

/**
 * \brief Registers the calculateAt method for a finite element class in Python.
 *
 * This function is used to expose the calculateAt method of a finite element class to Python using Pybind11.
 *
 * \tparam FE The finite element class type.
 * \tparam options Variadic template parameters for additional options when defining the Python class.
 * \param scope The Pybind11 handle to the Python module or class where the method will be registered.
 * \param cls The Pybind11 class wrapper for the finite element class.
 * \param restultTypesTuple A tuple containing the result types to be supported by the calculateAt method.
 *
 * \details
 * The calculateAt method is exposed to Python, allowing users to compute element values at a specific location.
 *
 * Example usage in Python:
 * \code{.py}
 * fe_instance = MyFiniteElementClass()
 * result = fe_instance.calculateAt(feRequirements, local, "linearStress")
 * \endcode
 *
 * \param feRequirements The requirements for the finite element calculation.
 * \param local The local coordinates where the calculation is performed.
 * \param resultType A string specifying the desired result type for the calculation.
 * \return The calculated result as a EigenType but NumPy array.
 *
 * \throws Dune::NotImplemented If the specified resultType is not supported by the finite element.
 */
template <class FE, class... options>
void registerCalculateAt(pybind11::handle scope, pybind11::class_<FE, options...> cls, auto resultTypesTuple) {
  using Traits         = typename FE::Traits;
  using FERequirements = typename FE::Requirement;
  cls.def(
      "calculateAt",
      [&](FE& self, const FERequirements& req, const Dune::FieldVector<double, Traits::mydim>& local,
          std::string resType) {
        Eigen::VectorXd result;
        bool success = false;
        Dune::Hybrid::forEach(resultTypesTuple, [&]<typename RT>(RT i) {
          if (resType == toString(i)) {
            success = true;
            result  = self.template calculateAt<RT::template Rebind>(req, local).asVec();
          }
        });
        if (success)
          return result;
        DUNE_THROW(Dune::NotImplemented, "Element " + Dune::className<FE>() + " doesn't support ResultType " + resType);
      },
      pybind11::arg("feRequirements"), pybind11::arg("local"), pybind11::arg("resultType"));
}

/**
 * \brief Register Python bindings for the FE class.
 *
 * This function registers Python bindings for a FE class, allowing it to be used in Python scripts.
 *
 * \tparam FE The FE class to be registered.
 * \tparam options Variadic template parameters for additional options when defining the Python class.
 *
 * \param scope A Pybind11 handle representing the Python scope where the class should be registered.
 * \param cls The Pybind11 class template to be used for registering the KirchhoffLoveShell class.
 *
 * \ingroup pythonbindings
 */
template <class FE, class... options>
void registerFE(pybind11::handle scope, pybind11::class_<FE, options...> cls) {
  using BH          = typename FE::BasisHandler;
  using GridElement = typename FE::GridElement;
  using Requirement = typename FE::Requirement;
  using FlatBasis   = typename FE::Traits::FlatBasis;

  int index = 0;
  cls.def(pybind11::init([](const BH& basisHandler, typename FE::PreTuple argsTuple) {
            auto unpackTuple = [&]<typename... Arg>(Arg&&... args) {
              return new FE(basisHandler, std::forward<Arg>(args)...);
            };
            return std::apply(unpackTuple, argsTuple);
          }),
          pybind11::keep_alive<1, 2>());

  cls.def("bind", [](FE& self, const GridElement& e) { self.bind(e); });
  cls.def("calculateScalar", [](FE& self, const Requirement& req, Ikarus::ScalarAffordance affordance) {
    return calculateScalar(self, req, affordance);
  });
  cls.def("calculateVector", [](FE& self, const Requirement& req, Ikarus::VectorAffordance affordance,
                                Eigen::Ref<Eigen::VectorXd> vec) { calculateVector(self, req, affordance, vec); });
  cls.def(
      "calculateMatrix",
      [](FE& self, const Requirement& req, Ikarus::MatrixAffordance affordance, Eigen::Ref<Eigen::MatrixXd> mat) {
        calculateMatrix(self, req, affordance, mat);
      },
      pybind11::arg("Requirement"), pybind11::arg("MatrixAffordance"), pybind11::arg("elementMatrix").noconvert());

  pybind11::module scopedf = pybind11::module::import("dune.functions");
  using LocalViewWrapper   = Dune::Python::LocalViewWrapper<FlatBasis>;

  auto includes    = Dune::Python::IncludeFiles{"dune/python/functions/globalbasis.hh"};
  auto [lv, isNew] = Dune::Python::insertClass<LocalViewWrapper>(
      scopedf, "LocalViewWrapper",
      Dune::Python::GenerateTypeName("Dune::Python::LocalViewWrapperWrapper", Dune::MetaType<FlatBasis>()), includes);

  if (isNew) {
    lv.def("bind", &LocalViewWrapper::bind);
    lv.def("unbind", &LocalViewWrapper::unbind);
    lv.def("index", [](const LocalViewWrapper& localView, int index) { return localView.index(index); });
    lv.def("__len__", [](LocalViewWrapper& self) -> int { return self.size(); });

    Dune::Python::Functions::registerTree<typename LocalViewWrapper::Tree>(lv);
    lv.def("tree", [](const LocalViewWrapper& view) { return view.tree(); });
  }

  cls.def(
      "localView",
      [](FE& self) -> LocalViewWrapper {
        auto lvWrapped = LocalViewWrapper(self.localView().globalBasis());
        // this can be simplified when https://gitlab.dune-project.org/staging/dune-functions/-/merge_requests/418
        // becomes available
        // pybind11::object obj = pybind11::cast(self.localView().element());
        lvWrapped.base().bind(self.localView().element());
        return lvWrapped;
      },
      pybind11::keep_alive<0, 1>());

  registerCalculateAt(scope, cls, typename FE::SupportedResultTypes());
  registerFERequirement(scope, cls);
}

} // namespace Ikarus::Python
