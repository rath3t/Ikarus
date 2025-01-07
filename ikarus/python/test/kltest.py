# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import debug_info

debug_info.setDebugFlags()
import ikarus as iks
from ikarus import finite_elements, assembler, io
import numpy as np
import scipy as sp

from dune.iga import ControlPoint, ControlPointNet, NurbsPatchData, IGAGrid
from dune.iga.basis import defaultGlobalBasis, Power, Nurbs

from dune.common.hashit import hashIt
from dune.iga.basis import preBasisTypeName
from dune.generator.generator import SimpleGenerator


def globalBasis(gv, tree):
    generator = SimpleGenerator("BasisHandler", "Ikarus::Python")

    pbfName = preBasisTypeName(tree, gv.cppTypeName)
    element_type = f"Ikarus::BasisHandler<{pbfName}>"
    includes = []
    includes += list(gv.cppIncludes)
    includes += ["dune/iga/nurbsbasis.hh"]
    includes += ["ikarus/python/basis/basis.hh"]

    moduleName = "Basis_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    basis = defaultGlobalBasis(gv, tree)
    return module.BasisHandler(basis)


if __name__ == "__main__":
    cp = ControlPoint((0, 0, 0), 1)
    cp2 = ControlPoint((10, 0, 0), 1)
    cp3 = ControlPoint((0, 2, 0), 1)
    cp4 = ControlPoint((10, 2, 0), 1)

    netC = ((cp, cp2), (cp3, cp4))
    net = ControlPointNet(netC)

    nurbsPatchData1 = NurbsPatchData(((0, 0, 1, 1), (0, 0, 1, 1)), net, (1, 1))
    nurbsPatchData2 = nurbsPatchData1.degreeElevate(0, 1)
    nurbsPatchDataFinal = nurbsPatchData2.degreeElevate(1, 1)
    # nurbsPatchDataFinal = nurbsPatchData2.knotRefinement((0.1,0.2,0.3,0.4,0.5), 0)
    # # nurbsPatchData = nurbsPatchData.knotRefinement((0.25,0.5,0.75), 1)

    gridView = IGAGrid(nurbsPatchDataFinal)

    gridView.hierarchicalGrid.globalRefine(2)
    basis = globalBasis(gridView, Power(Nurbs(), 3))

    flatBasis = basis.flat()
    print("Number of dofs: ", len(flatBasis))
    d = np.zeros(len(flatBasis))

    lambdaLoad = iks.Scalar(1.0)
    thickness = 0.1

    def vL(x, lambdaVal):
        return np.array([0, 0, 2 * thickness**3 * lambdaVal])

    vLoad = iks.finite_elements.volumeLoad3D(vL)

    klShell = iks.finite_elements.kirchhoffLoveShell(
        youngs_modulus=1000, nu=0.0, thickness=thickness
    )

    fes = []
    for e in gridView.elements:
        fes.append(iks.finite_elements.makeFE(basis, klShell, vLoad))
        fes[-1].bind(e)

    dirichletValues = iks.dirichletValues(flatBasis)

    def fixLeftAndRightEdge(vec, localIndex, localView, intersection):
        if (
            intersection.geometry.center[0] < 1e-8
            or intersection.geometry.center[0] > 10 - 1e-8
        ):
            vec[localView.index(localIndex)] = True

    dirichletValues.fixBoundaryDOFs(fixLeftAndRightEdge)

    assembler = iks.assembler.sparseFlatAssembler(fes, dirichletValues)

    def gradAndhess(dRedInput):
        req = fes[0].createRequirement()
        req.insertParameter(lambdaLoad)
        dBig = assembler.createFullVector(dRedInput)
        req.insertGlobalSolution(dBig)
        g = assembler.vector(req, iks.VectorAffordance.forces, iks.DBCOption.Reduced)
        h = assembler.matrix(req, iks.MatrixAffordance.stiffness, iks.DBCOption.Reduced)
        return [g, h]

    from numpy.linalg import norm

    maxiter = 100
    abs_tolerance = 1e-8
    d = np.zeros(assembler.reducedSize())
    for k in range(maxiter):
        R, K = gradAndhess(d)
        r_norm = norm(R)

        deltad = sp.sparse.linalg.spsolve(K, R)
        d -= deltad
        print(k, r_norm, norm(deltad))
        if r_norm < abs_tolerance:
            break

    dBig = assembler.createFullVector(d)
    displacementFunc = flatBasis.asFunction(dBig)
    vtkWriter = gridView.trimmedVtkWriter()
    vtkWriter.addPointData(displacementFunc, name="displacement")
    vtkWriter.write("KLshell")
    assert (
        0.2087577577980777 - max(d)
    ) < 1e-6, f"The maximum displacement should be 0.2087577577980777 but is {max(d)}"

    vtkWriter2 = iks.io.vtkWriter(assembler, io.DataCollector.iga)
    vtkWriter2.addPointData(displacementFunc, name="displacement")
    vtkWriter2.write("KLshell2")
