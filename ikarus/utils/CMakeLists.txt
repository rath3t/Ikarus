# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(
  FILES algorithms.hh
        autodiffhelper.hh
        basis.hh
        concepts.hh
        defaultfunctions.hh
        dirichletvalues.hh
        eigendunetransformations.hh
        eigensparseaddon.hh
        findlinesegment.hh
        flatprebasis.hh
        functionsanitychecks.hh
        functionhelper.hh
        init.hh
        linearalgebrahelper.hh
        makeenum.hh
        math.hh
        nonlinearoperator.hh
        nonlinopfactory.hh
        polyfit.hh
        pythonautodiffdefinitions.hh
        tensorutils.hh
        tensorproductquadrule.hh
        traits.hh
        traversal.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/utils
)

add_subdirectory(drawing)
add_subdirectory(observer)
