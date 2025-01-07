# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
add_subdirectory(ikarus)
set(CMAKE_CXX_STANDARD 20)

dune_python_configure_bindings(
  PATH
  "."
  PACKAGENAME
  ikarus
  CMAKE_METADATA_FLAGS
  DUNE_OPTS_FILE
  HAVE_MPI
  CMAKE_BUILD_TYPE
  CMAKE_CXX_COMPILER
  CMAKE_CXX_FLAGS
  CMAKE_C_COMPILER
  CMAKE_C_FLAGS
  CMAKE_PREFIX_PATH
  CMAKE_CXX_COMPILER_LAUNCHER
  CMAKE_CXX_STANDARD
)
if(POLICY CMP0087)
  cmake_policy(SET CMP0087 NEW)
endif()
