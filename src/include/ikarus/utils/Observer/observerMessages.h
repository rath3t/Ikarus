//
// Created by alex on 12/20/21.
//

#pragma once

enum class ControlMessages {
  BEGIN,
  CONTROL_STARTED,
  ITERATION_ENDED,
  LOADSTEP_ENDED,
  RESIDUALNORM_UPDATED,
  SOLUTION_CHANGED,
  END
};