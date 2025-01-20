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

struct Listener
{
  using Token = std::shared_ptr<void>;

  template <typename Broadcaster, typename F>
  void subscribe(Broadcaster& broadcaster, F&& f) {
    using Signature = typename traits::FunctionTraits<F>::FreeSignature;
    if constexpr (requires { broadcaster.operator->(); })
      t.push_back(broadcaster->template station<Signature>().registerListener(std::forward<F>(f)));
    else
      t.push_back(broadcaster.template station<Signature>().registerListener(std::forward<F>(f)));
  }

  template <typename Broadcaster, typename Signature, typename F>
  void subscribe(Broadcaster& broadcaster, F&& f) {
    if constexpr (requires { broadcaster.operator->(); })
      t.push_back(broadcaster->template station<Signature>().registerListener(std::forward<F>(f)));
    else
      t.push_back(broadcaster.template station<Signature>().registerListener(std::forward<F>(f)));
  }

  void unSubscribeAll() { t.clear(); }

private:
  std::vector<Token> t;
};

} // namespace Ikarus