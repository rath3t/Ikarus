// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/utils/makeenum.hh>

namespace Ikarus::Vtk {
/**
 * \enum DataTag
 * \brief Tag enum indicating cell data or point data.
 */
MAKE_ENUM(DataTag, asCellData, asPointData, asCellAndPointData);
} // namespace Ikarus::Vtk