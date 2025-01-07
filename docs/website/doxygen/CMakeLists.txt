# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-3.0-or-later
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS ../../../ikarus/**/*.hh ../../../ikarus/**/*.inl)
list(APPEND HEADERS modules.hh mainpage.md Doxylocal ../../literature.bib)
include(DuneDoxygen)
add_doxygen_target(DEPENDS ${HEADERS})
message(STATUS "Add doxygen target")
