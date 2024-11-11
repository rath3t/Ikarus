# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(
  FILES hyperelasticity.hh deviatoric.hh concepts.hh volumetric.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/finiteelements/mechanics/materials/hyperelasticity
)

add_subdirectory(modified)
add_subdirectory(regularized)
