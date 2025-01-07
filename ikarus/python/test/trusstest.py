# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later


import debug_info

debug_info.setDebugFlags()

import ikarus as iks
from ikarus import finite_elements, assembler
import numpy as np
import scipy as sp
import os

import dune.foamgrid
import dune.grid
import dune.functions


def trussTest(worldDim):
    assert worldDim == 2 or worldDim == 3
    gridDim = 1
    filedir = os.path.dirname(__file__)
    filename = os.path.join(filedir, f"auxiliaryfiles/truss{worldDim}d.msh")
    reader = (dune.grid.reader.gmsh, filename)

    grid = dune.foamgrid.foamGrid(reader, gridDim, worldDim)

    basis = iks.basis(
        grid, dune.functions.Power(dune.functions.Lagrange(order=1), worldDim)
    )

    flatBasis = basis.flat()
    d = np.zeros(len(flatBasis))

    lambdaLoad = iks.Scalar(0.0)

    fes = []

    tol = 1e-12

    if worldDim == 2:
        E = 30000
        A = 0.1
    else:
        E = 2000
        A = 1000

    trusses = finite_elements.truss(youngs_modulus=E, cross_section=A)

    for e in grid.elements:
        fes.append(finite_elements.makeFE(basis, trusses))
        fes[-1].bind(e)

    req = fes[0].createRequirement()
    req.insertParameter(lambdaLoad)

    assert req.parameter() == lambdaLoad
    req.insertGlobalSolution(d)

    d2 = req.globalSolution()

    # check that is really the same data address
    assert ("{}".format(hex(d2.__array_interface__["data"][0]))) == (
        "{}".format(hex(d.__array_interface__["data"][0]))
    )
    assert len(d2) == len(d)
    assert (d2 == d).all()

    forces = np.zeros(worldDim * 2)
    stiffness1 = np.zeros((worldDim * 2, worldDim * 2))
    stiffness2 = np.zeros((worldDim * 2, worldDim * 2))
    fes[0].calculateVector(req, iks.VectorAffordance.forces, forces)
    fes[0].calculateMatrix(req, iks.MatrixAffordance.stiffness, stiffness1)
    fes[1].calculateMatrix(req, iks.MatrixAffordance.stiffness, stiffness2)

    dirichletValues = iks.dirichletValues(flatBasis)

    def fixAllIndices(dirichletFlags, globalIndex):
        dirichletFlags[globalIndex] = True

    dirichletValues.fixBoundaryDOFs(fixAllIndices)
    if worldDim == 3:
        dirichletValues.setSingleDOF(1, True)

    denseAssembler = assembler.denseFlatAssembler(fes, dirichletValues)
    denseAssembler.bind(req, iks.AffordanceCollection.elastoStatics, iks.DBCOption.Full)

    Kdense = denseAssembler.matrix()
    forces = denseAssembler.vector()
    forces[4] += 1.0

    x = sp.linalg.solve(Kdense, -forces)

    if worldDim == 2:
        expectedDisp = (
            -0.169172906288868320
        )  # maximum vertical displacement of the center node
        expectedAxialForce = np.array([-4.592837713552849940, -4.592837713552849940])
    if worldDim == 3:
        expectedDisp = (
            -0.0008521190304919996
        )  # maximum vertical displacement of the center node
        expectedAxialForce = np.array([13.785353834067782586, -14.910220788494811472])

    assert (
        abs(abs(expectedDisp) - max(abs(x))) < tol
    ), f"The expected displacement should be {expectedDisp} but is {max(abs(x))}"

    req.insertGlobalSolution(x)

    # Test calculateAt Function
    N0 = fes[0].calculateAt(req, np.array([0.5]), "cauchyAxialForce")[0]
    N1 = fes[1].calculateAt(req, np.array([0.5]), "cauchyAxialForce")[0]
    N = np.array([N0, N1])

    for i in range(2):
        assert (
            abs(expectedAxialForce[i] - N[i]) < tol
        ), f"The expected Cauchy axial force for element {i} is {expectedAxialForce[i]} but is {N[i]}"

    # Querying for a different ResultType should result in a runtime error
    try:
        fes[0].calculateAt(req, np.array([0.5]), "PK2Stress")
    except RuntimeError:
        assert True
    else:
        assert False


if __name__ == "__main__":
    trussTest(2)
    trussTest(3)
