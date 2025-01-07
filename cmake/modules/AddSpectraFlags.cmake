# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# set HAVE_SPECTRA for config.h
set(HAVE_SPECTRA ${Spectra_FOUND})

# register all spectra related flags
if(Spectra_FOUND)
  dune_register_package_flags(LIBRARIES Spectra::Spectra COMPILE_DEFINITIONS "ENABLE_SPECTRA=1")
endif()

# add function to link against the spectra library
function(add_dune_Spectra_flags _targets)
  if(Spectra_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC Spectra::Spectra)
      target_compile_definitions(${_target} PUBLIC ENABLE_SPECTRA=1)
    endforeach(_target)
  endif()
endfunction(add_dune_Spectra_flags)
