# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import debug_info

debug_info.setDebugFlags()

import ikarus as iks
import dune.grid
import dune.functions

import math
import numpy as np


def makeGrid():
    lowerLeft = []
    upperRight = []
    elements = []
    for _ in range(2):
        lowerLeft.append(-1)
        upperRight.append(1)
        elements.append(3)

    return dune.grid.structuredGrid(lowerLeft, upperRight, elements)


def testDirichletValues():
    grid = makeGrid()
    basis = iks.basis(grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2))
    basis = basis.flat()

    dirichletValues = iks.dirichletValues(basis)

    assert basis.size() == dirichletValues.size
    assert basis.size() == len(dirichletValues)

    def fixOneIndex(vec, globalIndex):
        if globalIndex == 0:
            vec[globalIndex] = True

    dirichletValues.fixBoundaryDOFs(fixOneIndex)

    # This is equivalent to values[0] == True, but this syntax is discouraged by PEP
    assert dirichletValues.container[0]

    # Note that the result of localView.index(localIndex) is a multiIndex even for a flat basis, the localIndex is an int
    def fixAnotherIndexWithLocalView(vec, localIndex, localView):
        if localView.index(localIndex) == [4]:
            vec[localView.index(localIndex)] = True

        assert isinstance(localIndex, int)
        assert isinstance(localView.index(localIndex), list)
        assert isinstance(localView.index(localIndex)[0], int)

    dirichletValues.fixBoundaryDOFs(fixAnotherIndexWithLocalView)
    container = dirichletValues.container

    assert sum(container) == 2
    assert dirichletValues.fixedDOFsize == 2

    def fixTopSide(vec, localIndex, localView, intersection):
        if intersection.geometry.center[0] == 1.0:
            vec[localView.index(localIndex)] = True

        assert isinstance(localIndex, int)
        assert isinstance(localView.index(localIndex), list)
        assert isinstance(localView.index(localIndex)[0], int)

    dirichletValues.fixBoundaryDOFs(fixTopSide)

    # This assmues a structured grid
    indicesPerDirection: int = (math.sqrt(grid.size(0)) + 1) * 2
    assert dirichletValues.fixedDOFsize == 2 + indicesPerDirection

    # This checks whether container is a reference
    assert sum(container) == dirichletValues.fixedDOFsize

    # Test Subbasis
    dirichletValues2 = iks.dirichletValues(basis)

    dirichletValues2.fixBoundaryDOFs(fixOneIndex, 0)
    dirichletValues2.fixBoundaryDOFs(fixAnotherIndexWithLocalView, 0)
    dirichletValues2.fixBoundaryDOFs(fixTopSide, 0)

    assert dirichletValues2.fixedDOFsize == int(2 + indicesPerDirection / 2)

    dirichletValues2.fixBoundaryDOFs(fixTopSide, 1)
    assert dirichletValues2.fixedDOFsize == 2 + indicesPerDirection

    dirichletValues2.setSingleDOF(1, True)
    assert dirichletValues2.fixedDOFsize == 2 + indicesPerDirection + 1
    assert dirichletValues2.container[1]
    assert dirichletValues2.isConstrained(1)

    dirichletValues2.setSingleDOF((1), False)  # via MultiIndex
    assert dirichletValues2.fixedDOFsize == 2 + indicesPerDirection
    assert not dirichletValues2.isConstrained((1))  # via MultiIndex

    dirichletValues2.reset()
    assert dirichletValues2.fixedDOFsize == 0
    assert sum(dirichletValues2.container) == 0

    def fixDOFFunction(basis, vec):
        vec[:] = np.ones(dirichletValues2.size)

    dirichletValues2.fixDOFs(fixDOFFunction)
    assert dirichletValues2.fixedDOFsize == dirichletValues2.size


if __name__ == "__main__":
    testDirichletValues()
