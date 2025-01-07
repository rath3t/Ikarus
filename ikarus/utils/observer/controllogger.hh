// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controllogger.hh
 * \brief Observer implementation for logging control routines
 */

#pragma once
#include "observer.hh"
#include "observermessages.hh"

#include <chrono>

namespace Ikarus {
/**
 * \brief ControlLogger class for logging control messages.
 *
 * This class implements an observer for control messages and logs relevant information based on the received
 * messages.
 */
class ControlLogger : public IObserver<ControlMessages>
{
public:
  /**
   * \brief Implementation of the update method for control message logging.
   *
   * \param message The received control message.
   */
  void updateImpl(ControlMessages message) final;
  /**
   * \brief Implementation of the update method for logging control messages with string values.
   *
   * \param message The received control message.
   * \param val The string value associated with the message.
   */
  void updateImpl(ControlMessages message, const std::string& val) final;
  /**
   * \brief Implementation of the update method for logging control messages with an integer and a string value.
   *
   * \param message The received control message.
   * \param val1 The integer value associated with the message.
   * \param val2 The string value associated with the message.
   */
  void updateImpl(ControlMessages message, int val1, const std::string& val2) final;
  /**
   * \brief Implementation of the update method for logging control messages with an integer and a double value.
   *
   * \param message The received control message.
   * \param val1 The integer value associated with the message.
   * \param val2 The double value associated with the message.
   */
  void updateImpl(ControlMessages message, int val1, double val2) final;

private:
  using TimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
  TimePoint start_{};
  TimePoint stop_{};
  std::chrono::milliseconds duration_{};
};
} // namespace Ikarus
