# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.18)

file(
  GLOB_RECURSE sources
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  CONFIGURE_DEPENDS *.cpp
)

dune_add_library(
  ${PROJECT_NAME}
  ${sources}
  EXPORT_NAME
  Ikarus
  NAMESPACE
  "Ikarus::"
  LINK_LIBRARIES
  ${DUNE_LIBS}
)
target_link_dune_default_libraries(${PROJECT_NAME})

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)

add_subdirectory(assembler)
add_subdirectory(controlroutines)
add_subdirectory(finiteelements)
add_subdirectory(linearalgebra)
add_subdirectory(solver)
add_subdirectory(utils)
add_subdirectory(io)
if(DUNE_ENABLE_PYTHONBINDINGS)
  add_subdirectory(python)
endif()
