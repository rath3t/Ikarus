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

/**
 * \brief Implements a Broadcaster for a specific function signature with return type void.
 *
 * \tparam Args the arguments of the signature the Broadcaster can emit.
 */
template <typename... Args>
class Broadcaster<void(Args...)>
{
  using F = std::function<void(Args...)>;

  // The functions are stored as weak pointer, therefore the Broadcaster has no ownership over them.
  std::vector<std::weak_ptr<F>> listeners;

public:
  using Token = std::shared_ptr<F>;

  /**
   * \brief This method is used to register a Listener function.
   * \details The function that is passed in is first stored in a shared_ptr. After this, the shared_ptr is added to the
   * vector of listener functions, which leads to a implicit conversion to a weak_ptr. The shared_ptr is then returned
   * to the Listener that has called this function to be stored in a vector of shared_ptr<void> \ref listener.hh.
   *
   * \param target
   * \return Token
   */
  Token registerListener(F f) {
    auto sp = std::make_shared<F>(std::move(f));
    listeners.push_back(sp);
    return sp;
  }

  // Unregister a listener
  void unregisterListener(Token&& t) { t = nullptr; }

  /**
   * \brief This calls all the registered functions.
   */
  void notify(Args... args) {
    trim();
    for (auto& w : listeners) {
      if (auto p = w.lock()) {
        (*p)(args...);
      }
    }
  }

private:
  // Remove expired listeners, we need that because the weak pointers could already be invalidated (i.e. refcount == 0).
  void trim() {
    listeners.erase(std::remove_if(listeners.begin(), listeners.end(), [](auto& p) { return p.expired(); }),
                    listeners.end());
  }
};

/**
 * \brief Fuses together multiple function signatures that can be emitted by one broadcaster. A broadcaster has to be a
 * derived class of this class.
 *
 * \tparam Signatures
 */
template <typename... Signatures>
class Broadcasters : public Broadcaster<Signatures>...
{
public:
  using Broadcaster<Signatures>::registerListener...;
  using Broadcaster<Signatures>::unregisterListener...;
  using Broadcaster<Signatures>::notify...;

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

} // namespace Ikarus