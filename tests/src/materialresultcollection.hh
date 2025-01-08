// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <Eigen/Core>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/utils/makeenum.hh>

using namespace Ikarus;
using namespace Ikarus::Materials;

MAKE_ENUM(DeformationType, Undeformed, UniaxialTensile, BiaxialTensile, PureShear, Random);

template <typename ST, typename FD, typename SD>
auto initializeMaterialResults() {
  return std::make_tuple(ST{0.0}, FD::Zero().eval(), SD::Zero().eval());
}

template <DeformationType def>
auto OgdenTotalResults() {
  using DEV                                          = Ogden<3, PrincipalStretchTag::total>;
  using ST                                           = typename DEV::ScalarType;
  using FD                                           = typename DEV::FirstDerivative;
  using SD                                           = typename DEV::SecondDerivative;
  auto [energy, firstDerivatives, secondDerivatives] = initializeMaterialResults<ST, FD, SD>();

  if constexpr (def == DeformationType::UniaxialTensile) {
    energy = 30.11892638066439253852;
    firstDerivatives << -64.04041303584074143925, -64.04041303584074143925, 104.81337596088467916010;
    secondDerivatives.diagonal() << 585.26691789982063081054, 585.26691789982063081054, 136.38232022783703978279;
  } else if constexpr (def == DeformationType::BiaxialTensile) {
    energy = 102.74586311674921045255;
    firstDerivatives << -322.02260147539755251397, 104.81337596088467916010, 104.81337596088467916010;
    secondDerivatives.diagonal() << 1874.74758199925673974837, 136.38232022783703978279, 136.38232022783703978279;
  } else if constexpr (def == DeformationType::PureShear) {
    energy = 38.01387115054218563175;
    firstDerivatives << -137.80601871346158092286, 0.0, 104.81337596088467916010;
    secondDerivatives.diagonal() << 878.15470551965368326425, 379.33333333333333333333, 136.38232022783703978279;
  } else if constexpr (def == DeformationType::Random) {
    energy = 19.22492540140806116410;
    firstDerivatives << -149.31304210072199610872, -9.67967151876127425524, 0.0;
    secondDerivatives.diagonal() << 929.27792171584260115372, 407.61291331940502131823, 379.33333333333333333333;
  } else {
    secondDerivatives.diagonal() << 379.33333333333333333333, 379.33333333333333333333, 379.33333333333333333333;
  }

  return std::make_tuple(energy, firstDerivatives, secondDerivatives);
}

template <DeformationType def>
auto OgdenDeviatoricResults() {
  using DEV                                          = Ogden<3, PrincipalStretchTag::deviatoric>;
  using ST                                           = typename DEV::ScalarType;
  using FD                                           = typename DEV::FirstDerivative;
  using SD                                           = typename DEV::SecondDerivative;
  auto [energy, firstDerivatives, secondDerivatives] = initializeMaterialResults<ST, FD, SD>();

  if constexpr (def == DeformationType::UniaxialTensile) {
    energy = 30.11892638066439253852;
    firstDerivatives << -77.37108713810282996547, -77.37108713810282996547, 96.50011793037801997876;
    secondDerivatives << 505.66510597717411415639, -110.80811991756241388905, -133.28971931882038379099,
        -110.80811991756241388905, 505.66510597717411415639, -133.28971931882038379099, -133.28971931882038379099,
        -133.28971931882038379099, 25.36784593368447573792;
  } else if constexpr (def == DeformationType::BiaxialTensile) {
    energy = 102.74586311674921045255;
    firstDerivatives << -394.35652679503083661755, 76.68268938473846971177, 76.68268938473846971177;
    secondDerivatives << 2201.51629286968324379418, -140.23371403819243979711, -140.23371403819243979711,
        -140.23371403819243979711, 60.02021001397090733493, -117.42880594546643821184, -140.23371403819243979711,
        -117.42880594546643821184, 60.02021001397090733493;
  } else if constexpr (def == DeformationType::PureShear) {
    energy = 38.01387115054218563175;
    firstDerivatives << -157.44542092263587205377, -14.33533007968926359920, 94.34963137716988821177;
    secondDerivatives << 846.25615095796623372589, -113.99904353601931958853, -137.82141189493178989945,
        -113.99904353601931958853, 287.24803277351108881480, -128.00466564299765028820, -137.82141189493178989945,
        -128.00466564299765028820, 29.12771796706967717727;
  } else if constexpr (def == DeformationType::Random) {
    energy = 12.91548535350265217829;
    firstDerivatives << -108.96544359817898679151, 34.32155933680243633227, 44.29089744397990845154;
    secondDerivatives << 752.20966612104635712816, -160.54608411586555235210, -162.29649934275489849718,
        -160.54608411586555235210, 212.84045403016670728847, -161.61438739574503128460, -162.29649934275489849718,
        -161.61438739574503128460, 184.81860497880616173914;
  } else {
    secondDerivatives << 252.88888888888888888889, -126.44444444444444444444, -126.44444444444444444444,
        -126.44444444444444444444, 252.88888888888888888889, -126.44444444444444444444, -126.44444444444444444444,
        -126.44444444444444444444, 252.88888888888888888889;
  }

  return std::make_tuple(energy, firstDerivatives, secondDerivatives);
}