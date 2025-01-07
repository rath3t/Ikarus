# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator


def boundaryPatch(gridView, booleanVector):
    """
    @brief Creates a boundary patch for the given grid view and boolean vector.

    @param gridView: The grid view.
    @param booleanVector: The boolean vector.

    @return: The created boundary patch.
    """
    generator = MySimpleGenerator("BoundaryPatch", "Ikarus::Python")
    element_type = f"BoundaryPatch<{gridView.cppTypeName}>"

    includes = []
    includes += ["dune/fufem/boundarypatch.hh"]
    includes += ["ikarus/python/utils/boundarypatch.hh"]
    includes += gridView._includes
    moduleName = "boundaryPatch_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    return module.BoundaryPatch(gridView, booleanVector)


from io import StringIO
from dune.generator.algorithm import run
from dune.common import FieldVector


def globalIndexFromGlobalPosition(basis, pos):

    runCode = """
        #include <ikarus/utils/functionhelper.hh>
        #include <dune/python/pybind11/pybind11.h>
        #include <dune/python/pybind11/numpy.h>
        template <int size, typename Basis>
        auto callGlobalIndexFromGlobalPosition(const Basis& basis, const Dune::FieldVector<double, size>& pos) {
            auto arrayOfIndices = Ikarus::utils::globalIndexFromGlobalPosition(basis,pos);
            auto InnerSize = arrayOfIndices[0].size(); // How many indices are there? 1 for flat basis
            auto OuterSize = arrayOfIndices.size(); // How many dofs are there at the poss

            auto x = pybind11::array_t<size_t>({OuterSize,InnerSize});
            auto r = x.template mutable_unchecked<2>();
            for (pybind11::ssize_t i = 0; i < r.shape(0); i++)
                for (pybind11::ssize_t j = 0; j < r.shape(1); j++)
                    r(i,j) = arrayOfIndices[i][j];
            return x;
        }
    """
    pos = FieldVector(pos)

    return run("callGlobalIndexFromGlobalPosition", StringIO(runCode), basis, pos)
