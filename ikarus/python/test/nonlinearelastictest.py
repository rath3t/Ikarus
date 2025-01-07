# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import debug_info

debug_info.setDebugFlags()

import ikarus as iks
from ikarus import finite_elements, utils, assembler
import numpy as np
import scipy as sp
from scipy.optimize import minimize

import dune.grid
import dune.functions

if __name__ == "__main__":
    lowerLeft = []
    upperRight = []
    elements = []
    for i in range(2):
        lowerLeft.append(-1)
        upperRight.append(1)
        elements.append(3)

    grid = dune.grid.structuredGrid(lowerLeft, upperRight, elements)
    grid.hierarchicalGrid.globalRefine(0)
    basisLagrange1 = iks.basis(
        grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
    )
    flatBasis = basisLagrange1.flat()

    forces = np.zeros(8)
    stiffness = np.zeros((8, 8))

    def vL(x, lambdaVal):
        return np.array([lambdaVal * x[0] * 2, 2 * lambdaVal * x[1] * 0])

    vLoad = iks.finite_elements.volumeLoad2D(vL)

    def neumannLoad(x, lambdaVal):
        return np.array([lambdaVal * 0, lambdaVal])

    neumannVertices = np.zeros(grid.size(2), dtype=bool)

    def loadTopEdgePredicate(x):
        return True if x[1] > 0.9 else False

    indexSet = grid.indexSet
    for v in grid.vertices:
        neumannVertices[indexSet.index(v)] = loadTopEdgePredicate(v.geometry.center)

    boundaryPatch = iks.utils.boundaryPatch(grid, neumannVertices)

    nBLoad = iks.finite_elements.neumannBoundaryLoad(boundaryPatch, neumannLoad)

    svk = iks.materials.StVenantKirchhoff(E=1000, nu=0.3)

    svkPS = svk.asPlaneStress()

    nonLinElastic = iks.finite_elements.nonLinearElastic(svkPS)

    fes = []
    for e in grid.elements:
        fes.append(
            iks.finite_elements.makeFE(basisLagrange1, nonLinElastic, vLoad, nBLoad)
        )
        fes[-1].bind(e)

    dirichletValues = iks.dirichletValues(flatBasis)

    def fixFirstIndex(vec, globalIndex):
        vec[0] = True

    def fixAnotherVertex(vec, localIndex, localView):
        localView.index(localIndex)
        vec[1] = True

    def fixLeftHandEdge(vec, localIndex, localView, intersection):
        if intersection.geometry.center[1] < -0.9:
            vec[localView.index(localIndex)] = True

    dirichletValues.fixBoundaryDOFs(fixFirstIndex)
    dirichletValues.fixBoundaryDOFs(fixAnotherVertex)
    dirichletValues.fixBoundaryDOFs(fixLeftHandEdge)

    assembler = iks.assembler.sparseFlatAssembler(fes, dirichletValues)

    dRed = np.zeros(assembler.reducedSize())

    lambdaLoad = iks.Scalar(3.0)

    feReq = fes[0].createRequirement()

    def energy(dRedInput):
        feReq = fes[0].createRequirement()
        feReq.insertParameter(lambdaLoad)
        dBig = assembler.createFullVector(dRedInput)
        feReq.insertGlobalSolution(dBig)
        feReq.globalSolution()
        return assembler.scalar(feReq, iks.ScalarAffordance.mechanicalPotentialEnergy)

    def gradient(dRedInput):
        feReq = fes[0].createRequirement()
        feReq.insertParameter(lambdaLoad)
        dBig = assembler.createFullVector(dRedInput)
        feReq.insertGlobalSolution(dBig)
        return assembler.vector(
            feReq, iks.VectorAffordance.forces, iks.DBCOption.Reduced
        )

    def hess(dRedInput):
        feReq = fes[0].createRequirement()
        feReq.insertParameter(lambdaLoad)
        dBig = assembler.createFullVector(dRedInput)
        feReq.insertGlobalSolution(dBig)
        return assembler.matrix(
            feReq, iks.MatrixAffordance.stiffness, iks.DBCOption.Reduced
        ).todense()  # this is slow, but for this test we don't care

    resultd = minimize(energy, x0=dRed, options={"disp": True}, tol=1e-14)
    resultd2 = minimize(
        energy, x0=dRed, jac=gradient, options={"disp": True}, tol=1e-14
    )
    resultd3 = minimize(
        energy,
        method="trust-constr",
        x0=dRed,
        jac=gradient,
        hess=hess,
        options={"disp": True},
    )
    resultd4 = sp.optimize.root(gradient, jac=hess, x0=dRed, tol=1e-10)

    assert np.allclose(resultd.x, resultd2.x, atol=1e-6)
    assert np.allclose(resultd3.x, resultd4.x)
    assert np.all(abs(resultd3.grad) < 1e-8)
    assert np.all(abs(resultd4.fun) < 1e-8)

    feReq = fes[0].createRequirement()
    fullD = assembler.createFullVector(resultd2.x)
    feReq.insertGlobalSolution(fullD)

    res1 = fes[0].calculateAt(feReq, np.array([0.5, 0.5]), "PK2Stress")
