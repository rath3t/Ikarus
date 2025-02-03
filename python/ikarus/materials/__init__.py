# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt

from ikarus.generator import MySimpleGenerator

from enum import Enum
from typing import Union

# The list of supported dataCollectors
MuesliSmallStrain = Enum("MuesliSmallStrain", ["elasticIsotropicMaterial"])
MuesliFiniteStrain = Enum(
    "MuesliFiniteStrain",
    [
        "svkMaterial",
        "neohookeanMaterial",
        "mooneyMaterial",
        "yeohMaterial",
        "arrudaboyceMaterial",
    ],
)


def muesliMaterial(materialTag: Union[MuesliSmallStrain, MuesliFiniteStrain], **kwargs):
    assert isinstance(materialTag, MuesliSmallStrain) or isinstance(
        materialTag, MuesliFiniteStrain
    )

    materialTagSplitted = str(materialTag).split(".")
    wrappertype = (
        "SmallStrain" if isinstance(materialTag, MuesliSmallStrain) else "FiniteStrain"
    )
    materialtype = f"muesli::{materialTagSplitted[1]}"

    includes = []
    includes += ["ikarus/finiteelements/mechanics/materials/muesli/mueslimaterials.hh"]
    includes += ["ikarus/python/finiteelements/material.hh"]

    generator = MySimpleGenerator("MuesliMaterial", "Ikarus::Python")
    element_type = f"Ikarus::Materials::{wrappertype}<{materialtype}>"

    moduleName = "MuesliMaterial_" + hashIt(element_type)
    module = generator.load(
        includes=includes,
        typeName=element_type,
        moduleName=moduleName,
    )

    return module.MuesliMaterial(**kwargs)
