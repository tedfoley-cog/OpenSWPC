// Port of m_source.f90 (moment-tensor stress-glut and body-force sources).
#pragma once
#include <string>
#include <vector>

#include "readini.hpp"
#include "state.hpp"

namespace swpc {

struct Source {
    std::string stftype = "kupper";
    int nsrc = 0;
    int n_stfprm = 2;
    std::vector<std::vector<float>> srcprm;
    std::vector<float> sx, sz;
    std::vector<int> isrc, ksrc;
    std::vector<real_mp> mxx, mzz, mxz, fx, fz, mo;
    real_mp dt_dxz = 0;
    std::string sdep_fit = "asis";
};

void source_setup(State& s, Source& src, const IniFile& prm);
void source_stressglut(State& s, const Source& src, int it);
void source_bodyforce(State& s, const Source& src, int it);

}  // namespace swpc
