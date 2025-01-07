// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file genericobserver.hh
 * \brief Observer implementation for calling a callback function
 */

#pragma once
#include <concepts>
#include <string>

#include <spdlog/spdlog.h>

#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

/**
 * \brief GenericObserver class for observing specific messages.
 *
 * This class template implements an observer for a specific message type.
 *
 * \tparam M The type of messages to be observed.
 */
template <typename M>
class GenericObserver : public IObserver<M>
{
  using Messages = M;

public:
  /**
   * \brief Constructor for GenericObserver.
   *
   * Initializes the observer with a specific message and a function to be executed upon observation.
   *
   * \tparam F Type of the function to be executed.
   * \param p_message The message to be observed.
   * \param p_f The function to be executed with the current step .
   */
  template <typename F>
  GenericObserver(Messages message, F&& f)
      : message_{message},
        f_{f} {}
  void updateImpl(Messages message) override {
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
