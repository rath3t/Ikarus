# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
try:
    from dune.packagemetadata import registerExternalModule
    import pathlib

    # register ikarus to be recognized by dune-py (code generation module)
    # as a module of the dune universe
    registerExternalModule(
        moduleName="ikarus",
        modulePath=str(pathlib.Path(__file__).parent.resolve()),
    )

except ImportError:
    pass


import warnings
from ._ikarus import *
from .dirichlet_values import dirichletValues
from .basis import basis

from warnings import warn


def ValueWrapper(*args, **kwargs):
    warnings.warn(
        "The 'ValueWrapper' class was renamed to Scalar. Please use 'Scalar' instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return _ikarus.Scalar(*args, **kwargs)
