# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import types
from dune.common.hashit import hashIt
from .generator import MySimpleGenerator

from io import StringIO
from dune.generator.algorithm import run


def __fixBoundaryDOFs(dirichletValues):
    def __fixBoundaryDOFsFunc(dirichletValues, f, *args: int):
        prefixPathTypeName = "Dune::TypeTree::HybridTreePath<"
        prefixPathTypeName += ",".join(
            "Dune::index_constant<" + str(i) + ">" for i in args
        )
        prefixPathTypeName += ">"

        runCode = """
        #include <ikarus/python/dirichletvalues/dirichletvalues.hh> 
        template <typename DV>
        void fixBoundaryDofs(DV& dirichletValues, const pybind11::function& functor)
        {{
        Ikarus::Python::forwardCorrectFunction(dirichletValues, functor, 
        [&](auto&& functor_) {{ 
        return dirichletValues.fixBoundaryDOFs(functor_, {prefixPathType}{{}});
          }});
        }}
        """.format(
            prefixPathType=prefixPathTypeName
        )
        return run("fixBoundaryDofs", StringIO(runCode), dirichletValues, f)

    return __fixBoundaryDOFsFunc


def dirichletValues(basis):
    """
    @brief Creates a Dirichlet values handler for the given basis.

    @param basis: The basis.

    @return: The created Dirichlet values handler.
    """
    generator = MySimpleGenerator("DirichletValues", "Ikarus::Python")
    element_type = f"Ikarus::DirichletValues<{basis.cppTypeName},Eigen::VectorX<bool>>"

    includes = []
    includes += basis._includes
    includes += ["ikarus/utils/dirichletvalues.hh"]
    includes += ["ikarus/python/dirichletvalues/dirichletvalues.hh"]
    moduleName = "dirichletValues_" + hashIt(element_type)
    module = generator.load(
        includes=includes,
        typeName=element_type,
        moduleName=moduleName,
        dynamicAttr=True,
    )

    dv = module.DirichletValues(basis)
    dv.fixBoundaryDOFs = types.MethodType(__fixBoundaryDOFs(dv), dv)

    return dv
