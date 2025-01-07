# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator


def sparseFlatAssembler(fes, dirichletValues):
    """
    @brief Creates a sparse flat assembler for the given finite elements and Dirichlet values.

    @param fes: The list of finite elements.
    @param dirichletValues: The Dirichlet values.

    @return: The created sparse flat assembler.
    """
    element_type = f"Ikarus::SparseFlatAssembler<std::vector<{fes[0].cppTypeName}>,{dirichletValues.cppTypeName}>"
    generator = MySimpleGenerator("SparseFlatAssembler", "Ikarus::Python")

    includes = []
    includes += ["ikarus/assembler/simpleassemblers.hh"]
    includes += fes[0].cppIncludes  # include header of finite element
    includes += dirichletValues.cppIncludes
    includes += ["ikarus/python/assembler/flatassembler.hh"]
    moduleName = "SparseFlatAssembler_" + hashIt(element_type)
    module = generator.load(
        includes=includes,
        typeName=element_type,
        moduleName=moduleName,
        holder="std::shared_ptr",
    )
    return module.SparseFlatAssembler(fes, dirichletValues)


def denseFlatAssembler(fes, dirichletValues):
    """
    @brief Creates a dense flat assembler for the given finite elements and Dirichlet values.

    @param fes: The list of finite elements.
    @param dirichletValues: The Dirichlet values.

    @return: The created dense flat assembler.
    """
    element_type = f"Ikarus::DenseFlatAssembler<std::vector<{fes[0].cppTypeName}>,{dirichletValues.cppTypeName}>"
    generator = MySimpleGenerator("DenseFlatAssembler", "Ikarus::Python")

    includes = []
    includes += ["ikarus/assembler/simpleassemblers.hh"]
    includes += ["ikarus/python/assembler/flatassembler.hh"]
    includes += fes[0].cppIncludes  # include header of finite element
    includes += dirichletValues.cppIncludes
    moduleName = "DenseFlatAssembler_" + hashIt(element_type)
    module = generator.load(
        includes=includes,
        typeName=element_type,
        moduleName=moduleName,
        holder="std::shared_ptr",
    )
    return module.DenseFlatAssembler(fes, dirichletValues)


def assemblerManipulator(assembler):
    """
    @brief Creates an assembler manipulator for the given assembler. The manipulation happens by providing call back functions

    @param assembler: The assembler.

    @return: The created assembler manipulator.
    """
    element_type = f"decltype(Ikarus::makeAssemblerManipulator(std::declval<{assembler.cppTypeName}>()))::element_type"
    generator = MySimpleGenerator("AssemblerManipulator", "Ikarus::Python")

    includes = []
    includes += assembler.cppIncludes
    includes += ["ikarus/assembler/assemblermanipulatorfuser.hh"]
    includes += ["ikarus/python/assembler/flatassemblermanipulator.hh"]
    moduleName = "AssemblerManipulator_" + hashIt(element_type)
    module = generator.load(
        includes=includes,
        typeName=element_type,
        moduleName=moduleName,
        holder="std::shared_ptr",
    )
    return module.AssemblerManipulator(assembler)
