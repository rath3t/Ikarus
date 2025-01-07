# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(
  FILES interface.hh
        ogden.hh
  DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric
)
