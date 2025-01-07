# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# import debug_info
# debug_info.setDebugFlags()

import math

import ikarus as iks
from ikarus import finite_elements, assembler, io

import numpy as np
import scipy as sp

import dune.grid
import dune.functions

from dune.grid import gridFunction
from dune.vtk import FormatTypes, DataTypes

import os
import unittest
import pathlib

from dirichletvaluetest import makeGrid


class TestVtkWriter(unittest.TestCase):
    def setUp(self):
        reader = (
            dune.grid.reader.gmsh,
            os.path.join(os.path.dirname(__file__), "auxiliaryfiles/quad2d.msh"),
        )

        self.grid = dune.grid.ugGrid(reader, dimgrid=2)

        basis = iks.basis(
            self.grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
        )
        self.flatBasis = basis.flat()
        d = np.zeros(len(self.flatBasis))

        lambdaLoad = iks.Scalar(3.0)

        def vL(x, lambdaVal):
            return np.array([lambdaVal * x[0] * 2, 2 * lambdaVal * x[1] * 0])

        vLoad = finite_elements.volumeLoad2D(vL)

        fes = []
        linMat = iks.materials.LinearElasticity(E=1000, nu=0.2).asPlaneStress()
        linElastic = finite_elements.linearElastic(linMat)
        for e in self.grid.elements:
            fes.append(finite_elements.makeFE(basis, linElastic, vLoad))
            fes[-1].bind(e)

        req = fes[0].createRequirement()
        req.insertParameter(lambdaLoad)
        req.insertGlobalSolution(d)

        self.dirichletValues = iks.dirichletValues(self.flatBasis)

        def fixLeftHandEdge(vec, localIndex, localView, intersection):
            if intersection.geometry.center[0] < 1e-8:
                vec[localView.index(localIndex)] = True

        self.dirichletValues.fixBoundaryDOFs(fixLeftHandEdge)

        self.sparseAssembler = assembler.sparseFlatAssembler(fes, self.dirichletValues)
        self.sparseAssembler.bind(
            req, iks.AffordanceCollection.elastoStatics, iks.DBCOption.Full
        )

        Msparse = self.sparseAssembler.matrix()
        forces = self.sparseAssembler.vector()

        self.x = sp.sparse.linalg.spsolve(Msparse, -forces)
        req.insertGlobalSolution(self.x)

    def assertIsFile(self, path):
        if not pathlib.Path(path).resolve().is_file():
            raise AssertionError("File does not exist: %s" % str(path))

    def test(self):
        _ = io.vtkWriter(
            self.sparseAssembler,
            datatype=DataTypes.Float64,
            headertype=DataTypes.UInt32,
        )

        writer = io.vtkWriter(self.sparseAssembler, dataFormat=FormatTypes.ascii)
        writer.addInterpolation(
            self.x, self.flatBasis, "displacements", io.DataTag.asPointData
        )
        writer.addInterpolation(
            self.x, self.flatBasis, "displacements", io.DataTag.asCellData
        )

        xDisplacementBasis = dune.functions.subspaceBasis(self.flatBasis, 0)
        writer.addInterpolation(self.x, xDisplacementBasis, "u", io.DataTag.asPointData)

        writer.addAllResults(io.DataTag.asCellData)
        writer.addAllResults(io.DataTag.asPointData)

        writer.write("vtk_test")
        self.assertIsFile("vtk_test.vtu")

    def test_lagrange(self):
        writer2 = io.vtkWriter(
            self.sparseAssembler, dataCollector=io.DataCollector.lagrange, order=2
        )

        writer2.addResult("linearStress", io.DataTag.asCellData)
        writer2.addResult("linearStress", io.DataTag.asPointData)

        writer2.setFormat(FormatTypes.ascii)
        writer2.setDatatype(DataTypes.Float64)

        @gridFunction(self.grid)
        def g(x):
            return [math.sin(2 * math.pi * x[0] * x[1]), x[0] * x[1]] * 5

        writer2.addPointData(g, name="g", components=(0, 1))
        writer2.addPointData(g, name="g2", components=[0, 1, 2])

        writer2.write("vtk_test_lagrange")

    def test_structuredWriter(self):
        gridUG = makeGrid()

        basis = iks.basis(
            gridUG, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
        )
        flatBasisYASP = basis.flat()
        d = np.zeros(len(flatBasisYASP))

        lambdaLoad = iks.Scalar(3.0)

        fes = []
        linMat = iks.materials.LinearElasticity(E=1000, nu=0.2).asPlaneStress()
        linElastic = finite_elements.linearElastic(linMat)
        for e in gridUG.elements:
            fes.append(finite_elements.makeFE(basis, linElastic))
            fes[-1].bind(e)

        req = fes[0].createRequirement()
        lambdaLoad = iks.Scalar(3.0)
        req.insertParameter(lambdaLoad)
        req.insertGlobalSolution(d)

        dirichletValuesYASP = iks.dirichletValues(flatBasisYASP)
        sparseAssemblerYASP = assembler.sparseFlatAssembler(fes, dirichletValuesYASP)

        sparseAssemblerYASP.bind(
            req, iks.AffordanceCollection.elastoStatics, iks.DBCOption.Full
        )

        writer3 = io.vtkWriter(sparseAssemblerYASP, dataFormat=FormatTypes.ascii)
        writer3.addAllResults(io.DataTag.asCellData)

        writer3.write("vtk_test_structured")
        self.assertIsFile("vtk_test_structured.vtr")

        # self.assertEqual(fileName[-3:], "vtr")

    def test_discontinuousWriter(self):
        discontinuousVtkWriter = io.vtkWriter(
            self.sparseAssembler, dataCollector=io.DataCollector.discontinuous
        )
        discontinuousVtkWriter.addAllResults()  # defaults to pointData

        discontinuousVtkWriter.addInterpolation(
            self.x, self.flatBasis, "displacements", io.DataTag.asCellAndPointData
        )
        discontinuousVtkWriter.write("vtk_test_discontinous")


if __name__ == "__main__":
    unittest.main()
