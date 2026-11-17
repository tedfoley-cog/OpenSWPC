// Port of m_medium.f90, m_vmodel_uni.f90, m_vmodel_lhm.f90, m_seawater.f90.
#pragma once
#include "readini.hpp"
#include "state.hpp"

namespace swpc {

void medium_setup(State& s, const IniFile& prm);
float seawater_vel(float z_km, bool use_munk);

}  // namespace swpc
