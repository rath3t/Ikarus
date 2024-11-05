# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(
  FILES interface.hh
        linearelasticity.hh
        neohooke.hh
        blatzko.hh
        compressibleogden.hh
        incompressibleogden.hh
        deviatoric.hh
        volumetric.hh
        strainconversions.hh
        svk.hh
        tags.hh
        vanishingstress.hh
        vanishingstrain.hh
        vanishinghelpers.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/finiteelements/mechanics/materials
)
add_subdirectory(hyperelasticity)
