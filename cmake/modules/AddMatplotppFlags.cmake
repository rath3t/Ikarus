# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# set HAVE_MATPLOTLIB for config.h
set(HAVE_MATPLOTLIB ${Matplot++_FOUND})

# register all Matplot related flags
if(Matplot++_FOUND)
  dune_register_package_flags(
    LIBRARIES Matplot++::matplot COMPILE_DEFINITIONS "ENABLE_MATPLOTLIB=1"
  )
endif()

# add function to link against the Matplot library
function(add_dune_Matplot++_flags _targets)
  if(Matplot++_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC Matplot++::matplot)
      target_compile_definitions(${_target} PUBLIC ENABLE_MATPLOTLIB=1)
    endforeach(_target)
  endif()
endfunction(add_dune_Matplot++_flags)
