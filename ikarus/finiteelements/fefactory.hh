// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file fefactory.hh
 * \brief Contains the definition of the FEFactory class
 * \ingroup finiteelements
 */

#pragma once
#include <ikarus/finiteelements/febase.hh>
namespace Ikarus {
/**
 * \brief FEFactory is a convenient wrapper to forward arguments to PreFE and create and
 * construct a factory of finite elements.
 * \tparam BH The type of the basis handler.
 * \tparam SK Type of the skills.
 * \tparam useFlat A boolean indicating if the type of the underlying basis is of the flat or the untouched version.
 * \tparam useEigenRef A boolean flag indicating whether to use Eigen references.
 */
template <typename BH, typename SK, bool useFlat = true, bool useEigenRef = false>
struct FEFactory
{
  using Skills = SK;

private:
  const BH* basisHandler_;
  SK skills_;

public:
  /**
   * \brief constructor for FEFactory
   * \param basisHandler The basis handler.
   * \param sk Skill arguments.
   */
  template <typename SK2 = SK>
  FEFactory(const BH& basisHandler, const SK2& sk)
      : basisHandler_{&basisHandler},
        skills_{sk} {}

  auto operator()() {
    return std::apply(
        [&]<typename... Args>(Args&&... args) {
          // the template would not be needed,
          // when https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2019/p1814r0.html
          // will be implemented in clang. It is already implemented in gcc 12.2
          typename PreFE<BH, useFlat, useEigenRef>::template FE<std::decay_t<Args>::template Skill...> fe(
              *basisHandler_, std::forward<Args>(args)...);

          return fe;
        },
        skills_.args);
  }
};

/**
 * \brief A function to create a finite element using the flat version of the basis.
 * \tparam BH The type of the basis handler.
 * \tparam SK Type of the skills.
 * \tparam useFlat A boolean indicating if the  the underlying basis should be handed out as flat or the untouched
 * version.
 * \tparam useEigenRef A boolean flag indicating whether to use Eigen references.
 * \param basisHandler The basis handler.
 * \param sk Skill arguments.
 * \return An FEFactory object.
 */
template <bool useFlat = true, bool useEigenRef = false, typename BH, typename SK>
auto makeFE(const BH& basisHandler, const SK& sk) {
  FEFactory<BH, SK, useFlat, useEigenRef> factory(basisHandler, sk);

  return factory();
}

/**
 * \brief A function to create a finite element using the untouched version of the basis.
 * \tparam BH The type of the basis handler.
 * \tparam SK Type of the skills.

 * \tparam useEigenRef A boolean flag indicating whether to use Eigen references.
 * \param basisHandler The basis handler.
 * \param sk Skill arguments.
 * \return An FEFactory object.
 */
template <bool useEigenRef = false, typename BH, typename SK>
auto makeFEWithUnTouchedBasis(const BH& basisHandler, SK&& sk) {
  FEFactory<BH, SK, false, useEigenRef> factory(basisHandler, std::forward<SK>(sk));

  return factory();
}

} // namespace Ikarus
