# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# set HAVE_AUTODIFF for config.h
set(HAVE_AUTODIFF ${autodiff_FOUND})

# register all eigen related flags
if(autodiff_FOUND)
  dune_register_package_flags(LIBRARIES autodiff::autodiff COMPILE_DEFINITIONS "ENABLE_AUTODIFF=1")
endif()

# add function to link against the eigen library
function(add_dune_autodiff_flags _targets)
  if(autodiff_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC autodiff::autodiff)
      target_compile_definitions(${_target} PUBLIC ENABLE_AUTODIFF=1)
    endforeach(_target)
  endif()
endfunction(add_dune_autodiff_flags)
