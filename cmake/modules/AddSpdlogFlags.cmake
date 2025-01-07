# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# set HAVE_SPDLOG for config.h
set(HAVE_SPDLOG ${spdlog_FOUND})

# register all spdlog related flags
if(spdlog_FOUND)
  dune_register_package_flags(LIBRARIES spdlog::spdlog COMPILE_DEFINITIONS "ENABLE_SPDLOG=1")
endif()

# add function to link against the spdlog library
function(add_dune_spdlog_flags _targets)
  if(spdlog_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC spdlog::spdlog)
      target_compile_definitions(${_target} PUBLIC ENABLE_SPDLOG=1)
    endforeach(_target)
  endif()
endfunction(add_dune_spdlog_flags)
