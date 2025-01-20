# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(FILES controllogger.hh controlvtkwriter.hh genericlistener.hh nonlinearsolverlogger.hh
              listener.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/utils/listener
)
