// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file listener.hh
 * \brief Implementation of the observer design pattern with broadcasters
 */

#pragma once
#include <memory>
#include <vector>

#include <ikarus/utils/traits.hh>

namespace Ikarus {

/**
 * \brief
 * \details The functions that the listener is listening to are stored in a vector of shared_ptr<void>. This type
 * erasure has the advantage that we can listen to different function signatures and the Listener not beeing a template.
 * This works, because the deleter of the stored objects is not bound to the type information and thus there are no
 * memory leaks possible.
 *
 */
struct Listener
{
  using Token = std::shared_ptr<void>;

  /**
   * \brief Function to subscribe to a broadcaster with a given function (either a lambda, std::function or function
   * pointer).
   * \details This function deducts the types of the arguments itself and then forwards it. If there are problems with
   * type deduction use the function below and specify the types of the arguments manually.
   *
   * \tparam Broadcaster the type of the Broadcaster (for example a NonlinearSolver or ControlRoutine), can either be a
   * pointer or value.
   * \tparam F the type of the function
   * \param broadcaster the broadcaster
   * \param f the function
   */
  template <typename Broadcaster, typename F>
  void subscribe(Broadcaster& broadcaster, F&& f) {
    using Signature = typename traits::FunctionTraits<F>::FreeSignature;
    subscribe<Broadcaster, Signature>(broadcaster, std::forward<F>(f));
  }

  /**
   * \brief Function to subscribe to a broadcaster with a given function (either a lambda, std::function or function
   * pointer).
   * \tparam Broadcaster the type of the Broadcaster (for example a NonlinearSolver or ControlRoutine), can either be a
   * pointer or value.
   * \tparam Signature the exact signature of the function F
   * \tparam F the type of the function
   * \param broadcaster the broadcaster
   * \param f the function
   */
  template <typename Broadcaster, typename Signature, typename F>
  void subscribe(Broadcaster& broadcaster, F&& f) {
    if constexpr (requires { broadcaster.operator->(); })
      t.push_back(broadcaster->template station<Signature>().registerListener(std::forward<F>(f)));
    else
      t.push_back(broadcaster.template station<Signature>().registerListener(std::forward<F>(f)));
  }

  /**
   * \brief Unsubscribe from all listeners. At the moment unsubscibing can't be done more granularly.
   */
  void unSubscribeAll() { t.clear(); }

private:
  std::vector<Token> t;
};

} // namespace Ikarus