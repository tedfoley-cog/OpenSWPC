// Port of m_wav.f90 (station waveform sampling and SAC output).
#pragma once
#include <string>
#include <vector>

#include "readini.hpp"
#include "state.hpp"

namespace swpc {

struct Waveforms {
    int ntdec_w = 10, ntdec_w_prg = 0, ntw = 0, nst = 0;
    bool sw_v = false, sw_u = false, sw_stress = false, sw_strain = false;
    std::string wav_format = "sac";
    real_mp r40x = 0, r40z = 0, r41x = 0, r41z = 0;
    std::vector<float> xst, zst, stlo, stla;
    std::vector<std::string> stnm;
    std::vector<int> ist, kst;
    // wav_*[n][c][itw]
    std::vector<std::vector<std::vector<float>>> vel, disp, stress, strain;
    std::vector<float> ux, uz, exx, ezz, exz;
};

void wav_setup(State& s, Waveforms& w, const IniFile& prm);
void wav_store(State& s, Waveforms& w, int it);
void wav_write(const State& s, const Waveforms& w);

}  // namespace swpc
