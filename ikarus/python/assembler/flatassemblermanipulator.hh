// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file flatassemblerwrappers.hh
 * \brief Python bindings for assemblers
 */

#pragma once

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>
#include <dune/python/pybind11/stl_bind.h>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/python/assembler/flatassembler.hh>
#include <ikarus/python/finiteelements/scalarwrapper.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus::Python {

namespace Impl {

  template <typename NewType, typename F, typename... Args, size_t... Indices>
  decltype(auto) forward_last_as(F&& f, std::index_sequence<Indices...>, Args&&... args) {
    auto tup = std::forward_as_tuple(args...);
    return f(std::get<Indices>(tup)..., NewType(std::get<sizeof...(Args) - 1>(tup)));
  }

  template <typename NewType, typename F>
  auto wrapFunctionAndReplaceLastType(F&& f) {
    return [f](auto&&... args) -> decltype(auto) {
      return forward_last_as<NewType>(f, std::make_index_sequence<sizeof...(args) - 1>{},
                                      std::forward<decltype(args)>(args)...);
    };
  }
} // namespace Impl

// Since Pybind11 create a new scipy.sparse.csr_matrix from an Eigen::SparseMatrix, we have to create our own wrapper,
// which allows the modification of scalar entries of the sparse matrix in Python
template <typename T>
struct SparseMatrixWrapper
{
  SparseMatrixWrapper(Eigen::SparseMatrix<T>& matrix)
      : matrixRef(matrix) {}
  std::reference_wrapper<Eigen::SparseMatrix<T>> matrixRef;
};

template <typename T>
void registerSparseMatrixWrapper(pybind11::handle scope) {
  auto includes              = Dune::Python::IncludeFiles{"ikarus/finiteelements/flatassemblermanipulator.hh"};
  auto [lv, isNotRegistered] = Dune::Python::insertClass<SparseMatrixWrapper<T>>(
      scope, "SparseMatrixWrapper", Dune::Python::GenerateTypeName(Dune::className<SparseMatrixWrapper<T>>()),
      includes);
  if (isNotRegistered) {
    lv.def(pybind11::init<Eigen::SparseMatrix<T>&>())
        .def("__setitem__", [](SparseMatrixWrapper<T>& self, std::array<int, 2> index,
                               double val) { self.matrixRef.get().coeffRef(index[0], index[1]) = val; })
        .def("__getitem__", [](SparseMatrixWrapper<T>& self, std::array<int, 2> index) {
          return self.matrixRef.get().coeffRef(index[0], index[1]);
        });
  }
}

template <class AssemblerManipulator, class... options>
void registerAssemblerManipulator(pybind11::handle scope, pybind11::class_<AssemblerManipulator, options...> cls) {
  using pybind11::operator""_a;

  registerFlatAssembler(scope, cls);
  registerSparseMatrixWrapper<double>(scope);

  using UnderlyingAssembler = typename AssemblerManipulator::WrappedAssembler;

  cls.def(pybind11::init([](const UnderlyingAssembler& as) { return new AssemblerManipulator(as); }));

  using NewArgs = std::tuple<
      ScalarWrapper<std::reference_wrapper<typename AssemblerManipulator::ScalarType>>,
      Eigen::Ref<typename AssemblerManipulator::VectorType>,
      std::conditional_t<std::is_same_v<typename AssemblerManipulator::MatrixType, Eigen::SparseMatrix<double>>,
                         SparseMatrixWrapper<double>, Eigen::Ref<Eigen::MatrixXd>>>;
  Dune::Hybrid::forEach(Dune::Hybrid::integralRange(
                            Dune::index_constant<std::tuple_size_v<typename AssemblerManipulator::CallBackTypes>>()),
                        [&](auto i) {
                          using F          = std::tuple_element_t<i, typename AssemblerManipulator::CallBackTypes>;
                          using NewArg     = std::tuple_element_t<i, NewArgs>;
                          std::string name = std::string("add") +
                                             (i == 0     ? "Scalar"
                                              : (i == 1) ? "Vector"
                                                         : "Matrix") +
                                             std::string("CallBack");
                          using Traits            = Ikarus::traits::FunctionTraits<F>;
                          constexpr int lastIndex = Traits::numberOfArguments - 1;
                          // From Python we need a callback that accepts the wrapped types since otherwise Python
                          // creates copies and no modification is possible Therefore, from Python we get a callback in
                          // the style of Fmod= std::function<void(...,Wrapped<Type>)> but in the assembler we store can
                          // only store F=std::function<void(...,Type&)> wrapFunctionAndReplaceLastType takes care of
                          // this and wraps the "Fmod" call inside a "F" function
                          using FMod =
                              typename Ikarus::traits::ChangeArgTypeAtPos<F, lastIndex, NewArg>::NewFunctionType;

                          cls.def(name.c_str(), [&](AssemblerManipulator& self, FMod f) {
                            F fOrig = Impl::wrapFunctionAndReplaceLastType<NewArg>(std::forward<FMod>(f));

                            self.bind(std::move(fOrig));
                          });
                        });
}

} // namespace Ikarus::Python
