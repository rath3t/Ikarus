// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ikarus/utils/traits.hh>

struct A
{
  double foo(float b) const;
  double bar(int b, float c);
};

static int freeBar(double, float) { return 5; }

int main() {
  // Test lambda
  auto lambda = [](int a, double b) -> double { return a + b; };

  using TraitsLambda = Ikarus::traits::FunctionTraits<decltype(lambda)>;

  static_assert(std::is_same_v<typename TraitsLambda::args_type<0>, int>);
  static_assert(std::is_same_v<typename TraitsLambda::args_type<1>, double>);
  static_assert(std::is_same_v<typename TraitsLambda::return_type, double>);
  static_assert(TraitsLambda::numberOfArguments == 2);

  // Test free function
  freeBar(2.0, 7);
  using TraitsFree = Ikarus::traits::FunctionTraits<decltype(&freeBar)>;
  static_assert(std::is_same_v<typename TraitsFree::args_type<0>, double>);
  static_assert(std::is_same_v<typename TraitsFree::args_type<1>, float>);
  static_assert(std::is_same_v<typename TraitsFree::return_type, int>);
  static_assert(TraitsFree::numberOfArguments == 2);

  // Test const member function
  using TraitsConstMember = Ikarus::traits::FunctionTraits<decltype(&A::foo)>;
  static_assert(std::is_same_v<typename TraitsConstMember::args_type<0>, float>);
  static_assert(std::is_same_v<typename TraitsConstMember::return_type, double>);
  static_assert(TraitsConstMember::numberOfArguments == 1);

  // Test non-const member function
  using TraitsNonConstMember = Ikarus::traits::FunctionTraits<decltype(&A::bar)>;
  static_assert(std::is_same_v<typename TraitsNonConstMember::args_type<0>, int>);
  static_assert(std::is_same_v<typename TraitsNonConstMember::args_type<1>, float>);
  static_assert(std::is_same_v<typename TraitsNonConstMember::return_type, double>);
  static_assert(TraitsNonConstMember::numberOfArguments == 2);
}
