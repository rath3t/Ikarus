// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <iostream>

#include <ikarus/utils/init.hh>

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);

  std::cout << "Hello from ikarus!" << std::endl;
}
