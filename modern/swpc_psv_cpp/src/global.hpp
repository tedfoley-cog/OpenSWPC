// Port of m_global.f90 (single-process grid setup; MPI decomposition collapsed to nproc_x = 1).
#pragma once
#include "readini.hpp"
#include "state.hpp"

namespace swpc {
void global_setup(State& s, const IniFile& prm);
void global_setup2(State& s);
}  // namespace swpc
