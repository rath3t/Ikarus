# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(IkarusFormat)

if(ADD_FORMATTARGET)
  include(../CPM.cmake)
  CPMAddPackage(
    NAME Format.cmake
    VERSION 1.8.1
    GITHUB_REPOSITORY TheLartians/Format.cmake
    OPTIONS "FORMAT_SKIP_CMAKE NO" "CMAKE_FORMAT_EXCLUDE config.h.cmake"
  )
endif()

if(ADD_CODESPELLTARGET)
  find_package(Python)
  if(Python_FOUND)
    message("Adding codespell targets")
    set(CODESPELLCOMAND
        codespell
        -c
        -s
        -I
        codespellignore
        -S
        literature.bib
        --count
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/tests
        ${CMAKE_SOURCE_DIR}/examples
        ${CMAKE_SOURCE_DIR}/docs
        README.md
    )
    add_custom_target(
      codespell
      COMMAND ${CODESPELLCOMAND}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

    add_custom_target(
      fix-codespell
      COMMAND ${CODESPELLCOMAND} -w
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
  else()
    message(STATUS "python not found, codespell target not available")

  endif()

endif()

if(ADD_PYTHONFORMATTING)
  find_package(Python)
  if(Python_FOUND)
    message("Adding python format command")

    set(PYTHONFORMATTINGCOMAND
        pip
        install
        black==24.4.2
        &&
        python
        -m
        black
        .
    )
    add_custom_target(
      fix-pythonformat
      COMMAND ${PYTHONFORMATTINGCOMAND}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

    set(PYTHONFORMATTINGCOMANDCHECK
        pip
        install
        black==24.4.2
        &&
        python
        -m
        black
        .
        --check
    )
    add_custom_target(
      check-pythonformat
      COMMAND ${PYTHONFORMATTINGCOMANDCHECK}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
  else()
    message(STATUS "python not found, codespell target not available")
  endif()

endif()
