// Full-wavefield checkpoint dumps for state-level validation (same byte layout as the
// instrumented Fortran reference): enabled by SWPC_DUMP_DIR, cadence SWPC_DUMP_EVERY.
#pragma once
#include "state.hpp"

namespace swpc {
void dump_medium(const State& s);
void dump_state(const State& s, int it);
}  // namespace swpc
