# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# set HAVE_EIGEN for config.h
set(HAVE_EIGEN ${Eigen3_FOUND})

# register all eigen related flags
if(Eigen3_FOUND)
  dune_register_package_flags(LIBRARIES Eigen3::Eigen COMPILE_DEFINITIONS "ENABLE_EIGEN=1")
endif()

# add function to link against the eigen library
function(add_dune_Eigen3_flags _targets)
  if(Eigen3_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC Eigen3::Eigen)
      target_compile_definitions(${_target} PUBLIC ENABLE_EIGEN=1)
    endforeach(_target)
  endif()
endfunction(add_dune_Eigen3_flags)
