# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt

from ikarus.generator import MySimpleGenerator
from dune.vtk import FormatTypes, DataTypes
import dune.vtk

import types
from enum import Enum

# The list of supported dataCollectors
DataCollector = Enum("DataCollector", ["default", "lagrange", "discontinuous", "iga"])


def __addInterpolation(writer):
    def __addInterpolationFunc(
        writer, vals_, basis, name: str, dataTag=DataTag.asPointData
    ):
        gf = basis.asFunction(vals_)
        if dataTag == DataTag.asPointData or DataTag.asCellAndPointData:
            writer.addPointData(gf, name)
        if dataTag == DataTag.asCellData or DataTag.asCellAndPointData:
            writer.addCellData(gf, name)

    return __addInterpolationFunc


def vtkWriter(
    assembler,
    dataCollector: DataCollector = DataCollector.default,
    order: int = 1,
    dataFormat=FormatTypes.binary,
    datatype=DataTypes.Float32,
    headertype=DataTypes.UInt32,
):
    includes = []
    includes += ["dune/python/vtk/writer.hh"]
    includes += ["dune/vtk/writers/unstructuredgridwriter.hh"]
    includes += ["dune/vtk/writers/unstructuredgridwriter.hh"]
    includes += ["dune/vtk/datacollectors/lagrangedatacollector.hh"]
    includes += ["ikarus/io/vtkwriter.hh"]
    includes += ["ikarus/assembler/simpleassemblers.hh"]
    includes += assembler._includes

    gridViewName = assembler.gridView.cppTypeName
    dataCollectorName: str = ""

    if dataCollector is not DataCollector.default:
        if dataCollector == DataCollector.lagrange:
            dataCollectorName = (
                f"Dune::Vtk::LagrangeDataCollector<{gridViewName}, {order}>"
            )
        elif dataCollector == DataCollector.discontinuous:
            dataCollectorName = f"Dune::Vtk::DiscontinuousDataCollector<{gridViewName}>"
            includes += ["dune/vtk/datacollectors/discontinuousdatacollector.hh"]
        elif dataCollector == DataCollector.iga:
            dataCollectorName = (
                f"Dune::Vtk::DiscontinuousIgaDataCollector<{gridViewName}>"
            )
            includes += ["dune/iga/io/igadatacollector.hh"]

    generator = MySimpleGenerator("VtkWriter", "Ikarus::Python")

    defaultManager = f"Ikarus::Vtk::DefaultVTKWriterManager<{gridViewName}>"
    vtkWriterName = ""
    element_type = ""

    if dataCollector is DataCollector.default:
        vtkWriterName = f"typename {defaultManager}::template DefaultVTKWriter<>"
        dataCollectorName = f"typename {defaultManager}::DefaultDataCollector"
        element_type = f"Ikarus::Vtk::Writer<{assembler.cppTypeName}, {dataCollectorName}, {vtkWriterName}>"

    else:
        vtkWriterName = (
            f"typename {defaultManager}::template DefaultVTKWriter<{dataCollectorName}>"
        )
        element_type = f"Ikarus::Vtk::Writer<{assembler.cppTypeName}, {dataCollectorName}, {vtkWriterName}>"

    # Register VTKWriter
    dune.vtk.load(includes, vtkWriterName)

    includes += ["ikarus/python/io/vtkwriter.hh"]
    moduleName = "vtkWriter_" + hashIt(element_type)
    module = generator.load(
        includes=includes,
        typeName=element_type,
        moduleName=moduleName,
        baseClasses=[vtkWriterName],
        dynamicAttr=True,
    )
    writerModule = module.VtkWriter(assembler, dataFormat, datatype, headertype)
    writerModule.addInterpolation = types.MethodType(
        __addInterpolation(writerModule), writerModule
    )

    return writerModule
