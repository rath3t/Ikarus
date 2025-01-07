# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

option(DUNE_ENABLE_PYTHONBINDINGS "Enable Python bindings" ON)
option(ENABLE_WARNINGS "Enable warning" OFF)
option(ENABLE_TEST_COVERAGE "Gcov-type code coverage" OFF)
option(ENABLE_TEST_COVERAGE_SOURCEBASED "Source based code coverage" OFF)
option(DUNE_PYTHON_ALLOW_GET_PIP "Allow dune-common to install pip into venv" ON)

project(
  ikarus
  VERSION 0.4.1
  LANGUAGES C CXX
)

add_definitions(-DDUNE_LOCALFEFUNCTIONS_USE_EIGEN=1)
add_definitions(-DCMAKE_DISABLE_FIND_PACKAGE_Vc=1)
add_definitions(-DCMAKE_DISABLE_FIND_PACKAGE_LATEX=1)
add_definitions(-DBUILD_SHARED_LIBS=1)

# find dune packages to make cmake modules available
find_package(dune-common REQUIRED EXCLUDE_FROM_ALL)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules" ${dune-common_MODULE_PATH})
include(DuneMacros)

dune_project()
dune_enable_all_packages()

add_subdirectory(ikarus)
add_subdirectory(cmake)
if(NOT SKBUILD)
  add_subdirectory(sandbox)
  add_subdirectory(tests)
endif()

if(BUILD_DOCS)
  message("Build docs locally target ENABLED")
  add_subdirectory(docs)
endif()

if(DUNE_ENABLE_PYTHONBINDINGS)
  add_subdirectory(python)
endif()

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project()
