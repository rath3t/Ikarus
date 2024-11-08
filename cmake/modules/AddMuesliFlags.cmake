# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# set HAVE_MUESLI for config.h
set(HAVE_MUESLI true)

# register all muesli related flags
if(HAVE_MUESLI)
  message("Muesli found yuhu")
  dune_register_package_flags(LIBRARIES  muesli COMPILE_DEFINITIONS "ENABLE_MUESLI=1")
  else()
  message("Muesli not found crycry")
endif()



# add function to link against the muesli library
function(add_dune_muesli_flags _targets)
  if(HAVE_MUESLI)
    message("Muesli added")
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC  muesli)
      target_compile_definitions(${_target} PUBLIC ENABLE_MUESLI=1)

    endforeach(_target)
  endif()
endfunction(add_dune_muesli_flags)
