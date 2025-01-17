// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file genericlistener.hh
 * \brief Listener implementation for calling a callback function
 */

#pragma once
#include <spdlog/spdlog.h>

#include <ikarus/utils/broadcaster/listener.hh>
#include <ikarus/utils/observer/broadcastermessages.hh>

namespace Ikarus {

/**
 * \brief GenericListener class for observing specific messages. This class template implements an listener for a
 * specific message type.
 *
 * \tparam M The type of messages to be listend to.
 */
template <typename M>
class GenericListener : public Listener
{
  using Messages = M;

public:
  /**
   * \brief Constructor for GenericListener.
   *
   * Initializes the listener with a specific message and a function to be executed upon listening.
   *
   * \tparam F Type of the function to be executed.
   * \param message The message to be listend to.
   * \param f The function to be executed with the current step.
   */
  template <typename F>
  GenericListener(Messages message, F&& f)
      : message_{message},
        f_{f} {}

  template <typename BC>
  GenericListener& subscribeTo(BC& bc) {
    this->subscribe(bc, [&](M message) { this->updateImpl(message); });
    return *this;
  }

  void updateImpl(Messages message) {
    if (message_ == message) {
      f_(step_);
      ++step_;
    }
  }

private:
  Messages message_;
  std::function<void(int)> f_;
  int step_{0};
};
} // namespace Ikarus
