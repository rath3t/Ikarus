// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controllogger.hh
 * \brief Observer implementation for logging control routines
 */

#pragma once
#include <chrono>

#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/listener/listener.hh>

namespace Ikarus {
/**
 * \brief ControlLogger class for logging control messages.
 *
 * This class implements an observer for control messages and logs relevant information based on the received
 * messages.
 */
class ControlLogger : public Listener
{
public:
  template <typename BC>
  ControlLogger& subscribeTo(BC& bc) {
    this->subscribe(bc, [&](ControlMessages message) { this->updateImpl(message); });
    this->subscribe(bc, [&](ControlMessages message, const std::string& val) { this->updateImpl(message, val); });
    this->subscribe(
        bc, [&](ControlMessages message, int val1, const std::string& val2) { this->updateImpl(message, val1, val2); });
    this->subscribe(bc, [&](ControlMessages message, int val1, double val2) { this->updateImpl(message, val1, val2); });
    return *this;
  }

  /**
   * \brief Implementation of the update method for control message logging.
   *
   * \param message The received control message.
   */
  void updateImpl(ControlMessages message);
  /**
   * \brief Implementation of the update method for logging control messages with string values.
   *
   * \param message The received control message.
   * \param val The string value associated with the message.
   */
  void updateImpl(ControlMessages message, const std::string& val);
  /**
   * \brief Implementation of the update method for logging control messages with an integer and a string value.
   *
   * \param message The received control message.
   * \param val1 The integer value associated with the message.
   * \param val2 The string value associated with the message.
   */
  void updateImpl(ControlMessages message, int val1, const std::string& val2);
  /**
   * \brief Implementation of the update method for logging control messages with an integer and a double value.
   *
   * \param message The received control message.
   * \param val1 The integer value associated with the message.
   * \param val2 The double value associated with the message.
   */
  void updateImpl(ControlMessages message, int val1, double val2);

private:
  using TimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
  TimePoint start_{};
  TimePoint stop_{};
  std::chrono::milliseconds duration_{};
};
} // namespace Ikarus
