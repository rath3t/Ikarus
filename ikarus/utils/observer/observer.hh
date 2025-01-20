// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file observer.hh
 * \brief Implementation of the observer design pattern
 */

#pragma once
#include <map>
#include <memory>

#include <Eigen/Core>

#include <ikarus/utils/makeenum.hh>

namespace Ikarus {

/**
 * \brief Generic observer interface for the Observer design pattern.
 * See \cite gamma1995design for a description of the design pattern
 * \ingroup observer
 * \tparam MT The type of message that the observer can handle.
 */
template <typename MT>
class IObserver
{
public:
  using MessageType = MT;
  /** * \brief Virtual destructor for the observer interface. */
  virtual ~IObserver() = default;

  /**
   * \brief Update method for receiving notifications with a message.
   * \param message The message to be received.
   */
  void update(MessageType message) {
    checkMessageType(message);
    updateImpl(message);
  }

  /**
   * \brief Update method for receiving notifications with a message and a double value.
   * \param message The message to be received.
   * \param val The double value associated with the message.
   */
  void update(MessageType message, double val) {
    checkMessageType(message);
    updateImpl(message, val);
  }

  /**
   * \brief Update method for receiving notifications with a message and an integer value.
   * \param message The message to be received.
   * \param val The integer value associated with the message.
   */
  void update(MessageType message, int val) {
    checkMessageType(message);
    updateImpl(message, val);
  };

  /**
   * \brief Update method for receiving notifications with a message and a string value.
   * \param message The message to be received.
   * \param val The string value associated with the message.
   */
  void update(MessageType message, const std::string& val) {
    checkMessageType(message);
    updateImpl(message, val);
  };

  /**
   * \brief Update method for receiving notifications with a message and two values (integer and double).
   * \param message The message to be received.
   * \param val1 The integer value associated with the message.
   * \param val2 The double value associated with the message.
   */
  void update(MessageType message, int val1, double val2) {
    checkMessageType(message);
    updateImpl(message, val1, val2);
  };

  /**
   * \brief Update method for receiving notifications with a message, an integer value, and a string value.
   * \param message The message to be received.
   * \param val1 The integer value associated with the message.
   * \param val2 The string value associated with the message.
   */
  void update(MessageType message, int val1, const std::string& val2) {
    checkMessageType(message);
    updateImpl(message, val1, val2);
  };

  /**
   * \brief Update method for receiving notifications with a message and an Eigen::VectorXd.
   * \param message The message to be received.
   * \param vec The Eigen::VectorXd associated with the message.
   */
  void update(MessageType message, const Eigen::VectorXd& vec) {
    checkMessageType(message);
    updateImpl(message, vec);
  }

protected:
  virtual void updateImpl([[maybe_unused]] MessageType message) {};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] double val) {};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] int val) {};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] const std::string& val) {};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] int val1, const std::string& val2) {};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] int val1, double val2) {};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] const Eigen::VectorXd& vec) {}

private:
  void checkMessageType(MessageType message) {
    if (MessageType::END == message)
      DUNE_THROW(Dune::InvalidStateException, "The END enum type should not be used");
    if (MessageType::BEGIN == message)
      DUNE_THROW(Dune::InvalidStateException, "The BEGIN enum type should not be used");
  }
};

/**
 * \brief Generic observable interface for the Observer design pattern.
 * See \cite gamma1995design for a description of the design pattern
 * \tparam MessageType The type of message that the observable can handle.
 * \ingroup observer
 */
template <typename MessageType>
class IObservable
{
public:
  IObservable() {
    for (MessageType msg = MessageType::BEGIN; msg != MessageType::END; Ikarus::increment(msg))
      messages_.push_back(msg);
  }
  virtual ~IObservable() = default;
  /**
   * \brief Subscribe an observer to receive notifications for a specific message type.
   * \param message The message type to subscribe to.
   * \param observer The observer to be subscribed.
   */
  void subscribe(MessageType message, std::shared_ptr<IObserver<MessageType>> observer);
  /**
   * \brief Subscribe an observer to receive notifications for all message types.
   * \param observer The observer to be subscribed.
   */
  void subscribeAll(std::shared_ptr<IObserver<MessageType>> observer);
  /**
   * \brief Subscribe multiple observers to receive notifications for all message types.
   * \param observers List of observers to be subscribed.
   */
  void subscribeAll(std::initializer_list<std::shared_ptr<IObserver<MessageType>>> observers);
  /**
   * \brief Unsubscribe an observer from receiving notifications for a specific message type.
   * \param message The message type to unsubscribe from.
   * \param observer The observer to be unsubscribed.
   */
  void unSubscribe(MessageType message, std::shared_ptr<IObserver<MessageType>> observer);
  /**
   * \brief Unsubscribe an observer from receiving notifications for all message types.
   * \param observer The observer to be unsubscribed.
   */
  void unSubscribeAll(std::shared_ptr<IObserver<MessageType>> observer);
  /**
   * \brief Notify observers about a specific message type.
   * \param message The message type to notify about.
   */
  void notify(MessageType message);

  /**
   * \brief Notify observers about a specific message type with a floating-point value.
   * \tparam ScalarType The type of the floating-point value.
   * \param message The message type to notify about.
   * \param val The floating-point value to be sent with the notification.
   */
  template <std::floating_point ScalarType>
  void notify(MessageType message, ScalarType val);

  /**
   * \brief Notify observers about a specific message type with an integer value.
   * \param message The message type to notify about.
   * \param val The integer value to be sent with the notification.
   */
  void notify(MessageType message, int val);

  /**
   * \brief Notify observers about a specific message type with a string value.
   * \param message The message type to notify about.
   * \param val The string value to be sent with the notification.
   */
  void notify(MessageType message, const std::string& val);

  /**
   * \brief Notify observers about a specific message type with an integer and a double value.
   * \param message The message type to notify about.
   * \param val1 The integer value to be sent with the notification.
   * \param val2 The double value to be sent with the notification.
   */
  void notify(MessageType message, int val1, double val2);

  /**
   * \brief Notify observers about a specific message type with an integer value and a string value.
   * \param message The message type to notify about.
   * \param val1 The integer value to be sent with the notification.
   * \param val2 The string value to be sent with the notification.
   */
  void notify(MessageType message, int val1, const std::string& val2);

  /**
   * \brief Notify observers about a specific message type with an Eigen::VectorX.
   * \tparam ScalarType The type of the elements in the Eigen::VectorX.
   * \param message The message type to notify about.
   * \param vec The Eigen::VectorX to be sent with the notification.
   */
  template <std::floating_point ScalarType>
  void notify(MessageType message, Eigen::VectorX<ScalarType> vec);

private:
  using ObserverVector = std::vector<std::shared_ptr<IObserver<MessageType>>>;
  using ObserverMap    = std::map<MessageType, ObserverVector>;
  ObserverMap observers_;
  std::vector<MessageType> messages_;
};

template <typename MT>
void IObservable<MT>::subscribe(MT message, std::shared_ptr<IObserver<MT>> observer) {
  observers_[message];
  auto&& vectorOfObserversOfASpecificMessage = observers_[message];
  vectorOfObserversOfASpecificMessage.push_back(observer);
}

template <typename MT>
void IObservable<MT>::subscribeAll(std::shared_ptr<IObserver<MT>> observer) {
  for (auto& msg : messages_)
    subscribe(msg, observer);
}

template <typename MT>
void IObservable<MT>::subscribeAll(std::initializer_list<std::shared_ptr<IObserver<MT>>> observers) {
  for (auto& observer : observers)
    for (auto& msg : messages_)
      subscribe(msg, observer);
}

template <typename MT>
void IObservable<MT>::unSubscribe(MT message, std::shared_ptr<IObserver<MT>> observer) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  std::ranges::remove_if(vectorOfObserversOfASpecificMessage, [&observer](auto&& obs) { return obs == observer; });
}

template <typename MT>
void IObservable<MT>::unSubscribeAll(std::shared_ptr<IObserver<MT>> observer) {
  for (auto& msg : messages_)
    unSubscribe(msg, observer);
}

template <typename MT>
void IObservable<MT>::notify(MT message) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message);
}

template <typename MT>
template <std::floating_point ScalarType>
void IObservable<MT>::notify(MT message, ScalarType val) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, val);
}

template <typename MT>
void IObservable<MT>::notify(MT message, int val) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, val);
}

template <typename MT>
void IObservable<MT>::notify(MT message, const std::string& val) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, val);
}

template <typename MT>
void IObservable<MT>::notify(MT message, int val1, double val2) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, val1, val2);
}

template <typename MT>
void IObservable<MT>::notify(MT message, int val1, const std::string& val2) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, val1, val2);
}

template <typename MT>
template <std::floating_point ScalarType>
void IObservable<MT>::notify(MT message, Eigen::VectorX<ScalarType> vec) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, vec);
}
} // namespace Ikarus
