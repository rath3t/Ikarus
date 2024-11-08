# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

include(CMakeFindDependencyMacro)

find_package(spdlog)
include(AddSpdlogFlags)

find_package(Matplot++)
include(AddMatplotppFlags)

find_package(Eigen3 3.4.90 REQUIRED)
include(AddEigenFlags)

find_package(autodiff REQUIRED)
include(AddAutoDiffFlags)

find_package(Spectra REQUIRED)
include(AddSpectraFlags)

set(MUESLI_ROOT "/muesli")
add_library(muesli STATIC IMPORTED)
set_target_properties(muesli PROPERTIES
    IMPORTED_LOCATION "${MUESLI_ROOT}/lib/libmuesli_linux_ia64.a"
    INTERFACE_INCLUDE_DIRECTORIES "${MUESLI_ROOT}"
)

include(AddMuesliFlags)

