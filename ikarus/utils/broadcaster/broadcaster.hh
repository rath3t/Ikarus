// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file broadcaster.hh
 * \brief Implementation of the observer design pattern with broadcasters
 */

#pragma once
#include <functional>
#include <memory>
#include <vector>

namespace Ikarus {

template <class Args>
struct Broadcaster;

template <typename... Args>
class Broadcaster<void(Args...)>
{
  using f = std::function<void(Args...)>;
  std::vector<std::weak_ptr<f>> listeners;

public:
  using Token = std::shared_ptr<f>;

  // Register a listener
  Token registerListener(f target) {
    auto sp = std::make_shared<f>(std::move(target));
    listeners.push_back(sp);
    return sp;
  }

  // Unregister a listener
  void unregisterListener(Token&& t) { t = nullptr; }

  // Remove expired listeners
  void trim() {
    listeners.erase(std::remove_if(listeners.begin(), listeners.end(), [](auto& p) { return p.expired(); }),
                    listeners.end());
  }

  // Notify all listeners
  void notifyListeners(Args... args) {
    trim();
    for (auto& w : listeners) {
      if (auto p = w.lock()) {
        (*p)(args...);
      }
    }
  }
};

// BasicBroadcaster for multiple message types
template <typename... Ms>
class BasicBroadcaster : public Broadcaster<Ms>...
{
public:
  using Broadcaster<Ms>::registerListener...;
  using Broadcaster<Ms>::unregisterListener...;
  using Broadcaster<Ms>::notifyListeners...;

  // Access a specific broadcaster for a given message type
  template <typename M>
  Broadcaster<M>& station() {
    return *this;
  }

  // Access a specific broadcaster for a given message type (const version)
  template <typename M>
  const Broadcaster<M>& station() const {
    return *this;
  }
};

template <class... Signatures>
using Broadcasters = BasicBroadcaster<Signatures...>;

} // namespace Ikarus