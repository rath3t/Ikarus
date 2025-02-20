# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(
  FILES autodifffe.hh
        feconfiguration.hh
        fetraits.hh
        fehelper.hh
        physicshelper.hh
        febase.hh
        feresulttypes.hh
        mixin.hh
        fefactory.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/finiteelements
)

add_subdirectory(mechanics)
