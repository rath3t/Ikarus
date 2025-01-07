# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from dune.functions.globalbasis import preBasisTypeName
from dune.functions import defaultGlobalBasis
from .generator import MySimpleGenerator


def basis(gv, tree):
    """
    @brief Creates a basis handler for the given grid view and tree.

    @param gv: The grid view.
    @param tree: The tree of pre basis types.

    @return: The created basis handler.
    """
    generator = MySimpleGenerator("BasisHandler", "Ikarus::Python")

    pbfName = preBasisTypeName(tree, gv.cppTypeName)
    element_type = f"Ikarus::BasisHandler<{pbfName}>"

    includes = []
    includes += ["ikarus/python/basis/basis.hh"]
    includes += gv.cppIncludes
    moduleName = "Basis_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    basis = defaultGlobalBasis(gv, tree)
    return module.BasisHandler(basis)
