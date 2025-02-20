# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
set(IWYU_ARGS "-Xiwyu --mapping_file=iwyu.imp")
option(IKARUS_SANITIZERS "Add sanitizers to the compilation" OFF)
if(IKARUS_SANITIZERS)
  add_definitions(-DUSE_SANITIZER="Address,Memory,MemoryWithOrigins,Undefined,Thread,Leak")
endif(IKARUS_SANITIZERS)
include(../../cmake/tools.cmake)

file(
  GLOB programSourceFiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
   *.cpp
)
set(programSourceFiles
    testadaptivestepsizing.cpp
    testassembler.cpp
    testautodiffhelper.cpp
    testdependencies.cpp
    testdirichletvalue.cpp
    testenhancedassumedstrains.cpp
    testfebasis.cpp
    testferesulttypes.cpp
    testklshell.cpp
    testlinearelasticity.cpp
    testmanifolds.cpp
    testmaterial.cpp
    testnonlinearelasticityneohooke.cpp
    testnonlinearelasticitysvk.cpp
    testnonlinearoperator.cpp
    testpathfollowing.cpp
    testpolyfit.cpp
    testpythonconversion.cpp
    testshellthicknessquadrule.cpp
    testtrustregion.cpp
    testresultfunction.cpp
    testtruss.cpp
    testvtkwriter.cpp
    testhyperelasticity.cpp
    testmaterialswithad.cpp
)

set(TEST_DEPENDING_ON_LOCALFEFUNCTIONS
    testmanifolds
    testdirichletvalue
    testtrustregion
    testenhancedassumedstrains
    testassembler
    testlinearelasticity
    testnonlinearelasticityneohooke
    testnonlinearelasticitysvk
    testadaptivestepsizing
    testresultfunction.cpp
)
# Compile only tests
dune_add_test(
  SOURCES
  testfemixin.cpp
  LINK_LIBRARIES
  ikarus
  COMPILE_ONLY
  LABELS
  cpp
)
add_dune_pythonlibs_flags(testfemixin)
find_package(PythonLibs REQUIRED)
include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})
target_link_libraries(testfemixin PUBLIC ${PYTHON_LIBRARIES})
target_compile_features(testfemixin PUBLIC cxx_std_20)
dune_add_test(
  SOURCES
  testfunctiontraits.cpp
  LINK_LIBRARIES
  ikarus
  COMPILE_ONLY
  LABELS
  cpp
)
target_compile_features(testfunctiontraits PUBLIC cxx_std_20)

set(TEST_DEPENDING_ON_IGA testklshell testadaptivestepsizing)

set(TEST_NEED_MORE_TIME testadaptivestepsizing)

foreach(programSourceFile ${programSourceFiles})
  get_filename_component(programName ${programSourceFile} NAME_WLE)
  if(${programName} IN_LIST TEST_DEPENDING_ON_LOCALFEFUNCTIONS AND NOT HAVE_DUNE_LOCALFEFUNCTIONS)
    message(STATUS "Skipping test ${programName}, since dune-localfefunctions were not found")
    continue()
  endif()
  if(${programName} IN_LIST TEST_DEPENDING_ON_IGA AND NOT HAVE_DUNE_IGA)
    message(STATUS "Skipping test ${programName}, since dune-iga were not found")
    continue()
  endif()
  if(${programName} IN_LIST TEST_NEED_MORE_TIME)
    dune_add_test(SOURCES ${programSourceFile} LINK_LIBRARIES ikarus LABELS cpp)
    target_compile_features(${programName} PUBLIC cxx_std_20)
    set_tests_properties(${programName} PROPERTIES TIMEOUT 7500)
  else()
    dune_add_test(SOURCES ${programSourceFile} LINK_LIBRARIES ikarus LABELS cpp;cpp_quick)
    target_compile_features(${programName} PUBLIC cxx_std_20)
    set_tests_properties(${programName} PROPERTIES TIMEOUT 1000)
  endif()
  add_dune_pythonlibs_flags(${programName})
  find_package(PythonLibs REQUIRED)
  include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})
  target_link_libraries(${programName} PUBLIC ${PYTHON_LIBRARIES})

  target_compile_options(${programName} PUBLIC -ftemplate-backtrace-limit=0)
  if(ENABLE_WARNINGS)
    target_compile_options(
      ${programName}
      PUBLIC -Wall
             -Wextra
             -Wshadow
             -pedantic
             -Wnon-virtual-dtor
             -Wno-padded
             -Wno-c++20-compat
             -Wno-ctad-maybe-unsupported
             -Wno-sign-conversion
             -Wno-c++20-extensions
             -Wno-shorten-64-to-32
    )

    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      target_compile_options(
        ${programName}
        PUBLIC -Weverything -Wno-reserved-macro-identifier -Wno-exit-time-destructors
               -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-documentation
      )
    endif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  endif(ENABLE_WARNINGS)
  if(ENABLE_TEST_COVERAGE)
    target_compile_options(${programName} PUBLIC -O0 -g --coverage)
    target_link_options(${programName} PUBLIC --coverage)
    if(ENABLE_TEST_COVERAGE_SOURCEBASED)
      if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        target_compile_options(
          ${programName} PUBLIC -fno-limit-debug-info -fprofile-instr-generate -fcoverage-mapping
        )
        target_link_options(
          ${programName} PUBLIC -fno-limit-debug-info -fprofile-instr-generate -fcoverage-mapping
        )
      else()
        message(FATAL_ERROR "Source-based coverage is only available with clang")
      endif()
    endif()
  endif()

endforeach()

file(COPY testfiles/unstructuredtrianglesfine.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/testfiles)
