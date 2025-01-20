// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file listener.hh
 * \brief Implementation of the observer design pattern with broadcasters
 */

#pragma once
#include <memory>
#include <vector>

namespace Ikarus {

struct Listener
{
  using Token = std::shared_ptr<void>;

  // Control has to be derived from BroadCasters
  template <typename CONTROL, typename F>
  void subscribe(CONTROL& control, F&& f) {
    t.push_back(control.registerListener(f));
  }

  void unSubscribe() { t.clear(); }

private:
  std::vector<Token> t;
};

} // namespace Ikarus