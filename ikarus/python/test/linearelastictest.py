# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import debug_info

debug_info.setDebugFlags()

import ikarus as iks
from ikarus import finite_elements, utils, assembler
import numpy as np
import scipy as sp

import dune.grid
import dune.functions
from dune.vtk import vtkWriter, vtkUnstructuredGridWriter


def linElasticTest(easBool):
    lowerLeft = []
    upperRight = []
    elements = []
    for i in range(2):
        lowerLeft.append(-1)
        upperRight.append(1)
        elements.append(3)

    grid = dune.grid.structuredGrid(lowerLeft, upperRight, elements)
    # grid.hierarchicalGrid.globalRefine(4)
    basisLagrange12 = dune.functions.defaultGlobalBasis(
        grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
    )
    dirichletValuesT = iks.dirichletValues(basisLagrange12)

    basisLagrange1 = iks.basis(
        grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
    )
    flatBasis = basisLagrange1.flat()
    d = np.zeros(len(flatBasis))
    d[0] = 0.0

    lambdaLoad = iks.Scalar(3.0)

    fes = []

    def vL(x, lambdaVal):
        return np.array([lambdaVal * x[0] * 2, 2 * lambdaVal * x[1] * 0])

    vLoad = iks.finite_elements.volumeLoad2D(vL)

    def neumannLoad(x, lambdaVal):
        return np.array([lambdaVal * 0, lambdaVal])

    neumannVertices = np.zeros(grid.size(2), dtype=bool)

    def loadTopEdgePredicate(x):
        return True if x[1] > 0.999 else False

    indexSet = grid.indexSet
    for v in grid.vertices:
        neumannVertices[indexSet.index(v)] = loadTopEdgePredicate(v.geometry.center)

    boundaryPatch = iks.utils.boundaryPatch(grid, neumannVertices)

    nBLoad = iks.finite_elements.neumannBoundaryLoad(boundaryPatch, neumannLoad)

    linMat = iks.materials.LinearElasticity(E=1000, nu=0.2).asPlaneStress()

    linElastic = iks.finite_elements.linearElastic(linMat)
    easF = iks.finite_elements.eas(4)

    for e in grid.elements:
        if easBool:
            fes.append(
                iks.finite_elements.makeFE(
                    basisLagrange1, linElastic, easF, vLoad, nBLoad
                )
            )
        else:
            fes.append(
                iks.finite_elements.makeFE(basisLagrange1, linElastic, vLoad, nBLoad)
            )
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

    forces = np.zeros(8)
    stiffness = np.zeros((8, 8))
    fes[0].calculateVector(req, iks.VectorAffordance.forces, forces)
    fes[0].calculateMatrix(req, iks.MatrixAffordance.stiffness, stiffness)
    fes[0].localView()

    dirichletValues = iks.dirichletValues(flatBasis)

    def fixFirstIndex(vec, globalIndex):
        vec[0] = True

    def fixAnotherVertex(vec, localIndex, localView):
        localView.index(localIndex)
        vec[1] = True

    def fixLeftHandEdge(vec, localIndex, localView, intersection):
        if intersection.geometry.center[1] < -0.99999:
            vec[localView.index(localIndex)] = True

    dirichletValues.fixBoundaryDOFs(fixFirstIndex)
    dirichletValues.fixBoundaryDOFs(fixAnotherVertex)
    dirichletValues.fixBoundaryDOFs(fixLeftHandEdge)

    assembler = iks.assembler.sparseFlatAssembler(fes, dirichletValues)
    assemblerDense = iks.assembler.denseFlatAssembler(fes, dirichletValues)
    assembler.bind(req, iks.AffordanceCollection.elastoStatics, iks.DBCOption.Full)
    assemblerDense.bind(req, iks.AffordanceCollection.elastoStatics, iks.DBCOption.Full)

    Msparse = assembler.matrix()
    forces = assembler.vector()

    x = sp.sparse.linalg.spsolve(Msparse, -forces)
    fx = flatBasis.asFunction(x)
    # grid.plot()
    req.globalSolution()
    # Test calculateAt Function
    indexSet = grid.indexSet

    stressFuncScalar = grid.function(
        lambda e, x: fes[indexSet.index(e)].calculateAt(req, x, "linearStress")[0]
    )
    stressFuncVec = grid.function(
        lambda e, x: fes[indexSet.index(e)].calculateAt(req, x, "linearStress")[:]
    )
    # Writing results into vtk file

    fileName = "resultdisplacement" + ("EAS" if easBool else "")

    vtkWriter(grid, fileName, pointData={("displacement", (0, 1)): fx})

    writer2 = vtkUnstructuredGridWriter(grid)
    writer2.addCellData(stressFuncScalar, name="stress")
    writer2.addCellData(stressFuncVec, name="stress2")

    writer2.write(name="result" + ("EAS" if easBool else ""))

    # Querying for a different ResultType should result in a runtime error
    try:
        fes[0].calculateAt(req, np.array([0.5, 0.5]), "PK2Stress")
    except RuntimeError:
        assert True
    else:
        assert False

    assemblerManipulator = iks.assembler.assemblerManipulator(assembler)
    assemblerManipulatorDense = iks.assembler.assemblerManipulator(assemblerDense)

    if not easBool:

        def scalarf(assembler, req, affordance, scalar):
            scalar *= 2

        assemblerManipulator.addScalarCallBack(scalarf)
        assemblerManipulatorDense.addScalarCallBack(scalarf)
        req3 = assemblerManipulator.requirement()
        d3 = req3.globalSolution()
        d3 += 1
        assert (
            abs(2 * assembler.scalar() - assemblerManipulator.scalar()) < 1e-6
        ), f"assembler.scalar(): {assembler.scalar()}, assemblerManipulator.scalar(): {assemblerManipulator.scalar()}"
        assert (
            abs(2 * assembler.scalar() - assemblerManipulatorDense.scalar()) < 1e-6
        ), f"assembler.scalar(): {assembler.scalar()}, assemblerManipulator.scalar(): {assemblerManipulatorDense.scalar()}"

    def vectorf(assembler, req, affordance, dbcOption, vector):
        vector[5] += 2

    assemblerManipulator.addVectorCallBack(vectorf)
    assemblerManipulatorDense.addVectorCallBack(vectorf)
    assert (
        abs(assembler.vector()[5] - (assemblerManipulator.vector()[5] - 2)) < 1e-6
    ), f"assembler.vector()[5]: {assembler.vector()[5]}, assemblerManipulator.vector()[5]: {assemblerManipulator.vector()[5]}"
    assert (
        abs(assembler.vector()[5] - (assemblerManipulatorDense.vector()[5] - 2)) < 1e-6
    ), f"assembler.vector()[5]: {assembler.vector()[5]}, assemblerManipulator.vector()[5]: {assemblerManipulatorDense.vector()[5]}"

    def matrixf(assembler, req, affordance, dbcOption, matrix):
        matrix[5, 6] = matrix[5, 6] + 2.0

    assemblerManipulator.addMatrixCallBack(matrixf)
    assemblerManipulatorDense.addMatrixCallBack(matrixf)
    assert (
        abs(assembler.matrix()[5, 6] - (assemblerManipulator.matrix()[5, 6] - 2)) < 1e-6
    ), f"assembler.matrix()[5,6]: {assembler.matrix()[5,6]}, assemblerManipulator.matrix()[5,6]: {assemblerManipulator.matrix()[5,6]}"
    assert (
        abs(assembler.matrix()[5, 6] - (assemblerManipulatorDense.matrix()[5, 6] - 2))
        < 1e-6
    ), f"assembler.matrix()[5,6]: {assembler.matrix()[5,6]}, assemblerManipulator.matrix()[5,6]: {assemblerManipulatorDense.matrix()[5,6]}"

    startX = -1
    incX = 2 / 3
    startY = -1
    incY = 2 / 3
    for j in range(4):
        for i in range(4):
            indices = iks.utils.globalIndexFromGlobalPosition(
                flatBasis, (startX + i * incX, startY + j * incY)
            )
            assert (
                indices[0] == i + 4 * j
            ), f"Expected {i+4*j} from {i},{j}, got {indices[0]} at pos {(startX+i*incX,startY+j*incY)}"
            assert (
                indices[1] == i + 4 * j + 16
            ), f"Expected {i+4*j+16} from {i},{j}, got {indices[1]} at pos {(startX+i*incX,startY+j*incY)}"

    assert assembler.size == assemblerManipulator.size == flatBasis.size()
    assert len(assembler) == len(assemblerManipulator) == len(flatBasis)

    assert (
        assembler.isConstrained(1)
        == assemblerManipulator.isConstrained(1)
        == dirichletValues.isConstrained(1)
    )
    assert assembler.constraintsBelow(i) == assemblerManipulator.constraintsBelow(i)


if __name__ == "__main__":
    linElasticTest(easBool=False)
    linElasticTest(easBool=True)
