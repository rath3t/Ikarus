# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import debug_info

debug_info.setDebugFlags()

import ikarus as iks
from ikarus import finite_elements, utils

import numpy as np

import dune.grid
import dune.functions

import unittest

from dirichletvaluetest import makeGrid


class DynamicsTest(unittest.TestCase):
    def setUp(self):
        self.grid = makeGrid()

        basis = iks.basis(
            self.grid,
            dune.functions.Power(
                dune.functions.Lagrange(order=1), 2, layout="interleaved"
            ),
        )
        self.flatBasis = basis.flat()

        self.fes = []
        linMat = iks.materials.LinearElasticity(E=1000, nu=0.2).asPlaneStrain()
        linElastic = finite_elements.linearElastic(linMat, density=100)
        for e in self.grid.elements:
            self.fes.append(finite_elements.makeFE(basis, linElastic))
            self.fes[-1].bind(e)

        self.dirichletValues = iks.dirichletValues(self.flatBasis)

        def fixLeftHandEdge(vec, localIndex, localView, intersection):
            if intersection.geometry.center[0] < 1e-8:
                vec[localView.index(localIndex)] = True

        self.dirichletValues.fixBoundaryDOFs(fixLeftHandEdge)

    def test(self):
        mA = utils.modalAnalysis(self.fes, self.dirichletValues)
        mA.compute()
        frequencies = mA.angularFrequencies()

        mA.bindLumpingScheme()
        mA.compute()
        frequenciesLumped = mA.angularFrequencies()

        self.assertTrue(np.all(frequencies >= frequenciesLumped))

        mA.unBindLumpingScheme()
        mA.compute()
        frequenciesNotLumpedAnymore = mA.angularFrequencies()

        self.assertTrue(np.isclose(frequencies, frequenciesNotLumpedAnymore).all())


if __name__ == "__main__":
    unittest.main()
