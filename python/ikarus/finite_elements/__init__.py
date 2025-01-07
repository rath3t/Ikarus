# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator


def registerPreElement(name, includes, element_type, *args):
    """
    @brief Registers a pre-element with the specified parameters.

    @param name: The name of the pre-element.
    @param includes: List of additional include files.
    @param element_type: The type of the pre-element.
    @param args: Additional arguments required for the pre-element registration.

    @return: The registered pre-element function.
    """
    generator = MySimpleGenerator(name, "Ikarus::Python")
    includes += ["ikarus/python/finiteelements/registerpreelement.hh"]

    moduleName = name + "_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    f = getattr(module, name)
    return f(*args)


def volumeLoad(f, d: int):
    """
    @brief Creates a volume load pre-element for the specified dimension.

    @param f: The volume load function.
    @param d: The dimension of the volume load (1, 2, or 3).

    @return: The registered volume load pre-element function.
    """
    includes = ["ikarus/finiteelements/mechanics/loads/volume.hh"]
    includes += ["dune/python/pybind11/functional.h"]
    includes += ["dune/python/pybind11/stl.h"]
    includes += ["dune/python/pybind11/eigen.h"]

    element_type = f"Ikarus::VolumeLoadPre<{d}>"
    return registerPreElement("VolumeLoadPre", includes, element_type, f)


def volumeLoad1D(f):
    """
    @brief Creates a 1D volume load pre-element.

    @param f: The volume load function.

    @return: The registered 1D volume load pre-element function.
    """
    return volumeLoad(f, 1)


def volumeLoad2D(f):
    """
    @brief Creates a 2D volume load pre-element.

    @param f: The volume load function.

    @return: The registered 2D volume load pre-element function.
    """
    return volumeLoad(f, 2)


def volumeLoad3D(f):
    """
    @brief Creates a 3D volume load pre-element.

    @param f: The volume load function.

    @return: The registered 3D volume load pre-element function.
    """
    return volumeLoad(f, 3)


def neumannBoundaryLoad(boundaryPatch, f):
    """
    @brief Creates a Neumann boundary load pre-element.

    @param boundaryPatch: The boundary patch.
    @param f: The traction function.

    @return: The registered Neumann boundary load pre-element function.
    """
    includes = ["ikarus/finiteelements/mechanics/loads/traction.hh"]
    includes += ["dune/python/pybind11/functional.h"]
    includes += ["dune/python/pybind11/stl.h"]
    includes += ["dune/python/pybind11/eigen.h"]
    includes += boundaryPatch._includes
    element_type = (
        f"Ikarus::NeumannBoundaryLoadPre<{boundaryPatch.gridView().cppTypeName}>"
    )
    return registerPreElement(
        "NeumannBoundaryLoadPre", includes, element_type, boundaryPatch, f
    )


def nonLinearElastic(mat):
    """
    @brief Creates a non-linear elastic pre-element.

    @param mat: The material.

    @return: The registered non-linear elastic pre-element function.
    """
    includes = ["ikarus/finiteelements/mechanics/nonlinearelastic.hh"]

    element_type = f"Ikarus::NonLinearElasticPre<{mat.cppTypeName}>"
    includes += mat._includes
    return registerPreElement("NonLinearElasticPre", includes, element_type, mat)


def linearElastic(mat):
    """
    @brief Creates a linear elastic pre-element.

    @param mat: The material.

    @return: The registered linear elastic pre-element function.
    """
    includes = ["ikarus/finiteelements/mechanics/linearelastic.hh"]
    element_type = f"Ikarus::LinearElasticPre<{mat.cppTypeName}>"
    includes += mat._includes
    return registerPreElement("LinearElasticPre", includes, element_type, mat)


def truss(youngs_modulus, cross_section):
    """
    @brief Creates a truss pre-element.

    @param youngs_modulus: Young's modulus.
    @param cross_section: Cross-section area.

    @return: The registered truss pre-element function.
    """
    includes = ["ikarus/finiteelements/mechanics/truss.hh"]
    element_type = "Ikarus::TrussPre"
    return registerPreElement(
        "TrussPre", includes, element_type, youngs_modulus, cross_section
    )


def eas(numberofparameters):
    """
    @brief Creates an enhanced assumed strains pre-element.

    @param numberofparameters: The number of parameters.

    @return: The registered enhanced assumed strains pre-element function.
    """
    includes = ["ikarus/finiteelements/mechanics/enhancedassumedstrains.hh"]
    element_type = "Ikarus::EnhancedAssumedStrainsPre"
    return registerPreElement(
        "EnhancedAssumedStrainsPre", includes, element_type, numberofparameters
    )


def kirchhoffLoveShell(youngs_modulus: float, nu, thickness):
    """
    @brief Creates a Kirchhoff Love shell pre-element.

    @param youngs_modulus: Young's modulus.
    @param nu: Poisson's ratio.
    @param thickness: The thickness of the shell.

    @return: The registered Kirchhoff Love shell pre-element function.
    """
    includes = ["ikarus/finiteelements/mechanics/kirchhoffloveshell.hh"]
    element_type = "Ikarus::KirchhoffLoveShellPre"
    return registerPreElement(
        "KirchhoffLoveShellPre", includes, element_type, youngs_modulus, nu, thickness
    )


def makeFE(basis, *skills):
    """
    @brief Creates a finite element with the specified basis and skills.

    @param basis: The basis for the finite element.
    @param skills: Skills for the finite element.

    @return: The created finite element.
    """
    includes = ["ikarus/python/finiteelements/fe.hh"]
    preFE_type = f"Ikarus::PreFE<{basis.cppTypeName},true,true>"
    element_type = f"Ikarus::FE<{preFE_type},"

    for arg in skills:
        element_type += arg.cppTypeName + "::Skill,"
    element_type = element_type[:-1]
    element_type += ">"

    generator = MySimpleGenerator("FE", "Ikarus::Python")

    includes += basis._includes
    for arg in skills:
        includes += arg._includes

    moduleName = "FE" + "_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    tupleskills = tuple(skills)
    return module.FE(basis, tupleskills)
