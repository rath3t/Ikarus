# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

add_subdirectory(assembler)
add_subdirectory(finite_elements)
add_subdirectory(utils)
add_subdirectory(io)

add_python_targets(ikarus __init__ dirichlet_values basis generator)

dune_add_pybind11_module(NAME _ikarus)

set_property(
  TARGET _ikarus
  PROPERTY LINK_LIBRARIES dunecommon ikarus
  APPEND
)

if(SKBUILD)
  install(TARGETS _ikarus LIBRARY DESTINATION python/ikarus)
endif()
