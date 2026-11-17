#include <cmath>
#include <ctime>

#include "fdtool.hpp"
#include "sac.hpp"
#include "wav.hpp"

namespace swpc {

namespace {
SacHeader base_header(const State& s, const Waveforms& w, int n, float mag) {
    SacHeader h;
    h.f[EVLO - 1] = s.evlo;
    h.f[EVLA - 1] = s.evla;
    h.f[EVDP - 1] = s.evdp;
    h.f[B - 1] = s.tbeg;
    h.f[DELTA - 1] = w.ntdec_w * s.dt;
    h.f[MAG - 1] = mag;
    h.i[NPTS - 71] = w.ntw;
    if (s.bf_mode) {
        h.f[USER0 - 1] = s.fx0;
        h.f[USER2 - 1] = s.fz0;
    } else {
        h.f[USER0 - 1] = s.mxx0;
        h.f[USER2 - 1] = s.mzz0;
        h.f[USER4 - 1] = s.mxz0;
    }
    h.f[USER6 - 1] = s.clon;
    h.f[USER7 - 1] = s.clat;
    h.f[USER8 - 1] = s.phi;
    h.f[O - 1] = s.otim;
    std::time_t t = s.exedate;
    std::tm lt{};
    localtime_r(&t, &lt);
    h.i[NZYEAR - 71] = lt.tm_year + 1900;
    h.i[NZJDAY - 71] = lt.tm_yday + 1;
    h.i[NZHOUR - 71] = lt.tm_hour;
    h.i[NZMIN - 71] = lt.tm_min;
    h.i[NZSEC - 71] = lt.tm_sec;
    h.i[NZMSEC - 71] = 0;
    h.k[KEVNM] = s.title.substr(0, 16);
    h.k[KSTNM] = w.stnm[n];
    h.f[STLO - 1] = w.stlo[n];
    h.f[STLA - 1] = w.stla[n];
    h.f[STDP - 1] = w.zst[n] * 1000;
    h.l[LCALDA - 106] = 0;
    h.f[DIST - 1] = std::sqrt((s.sx0 - w.xst[n]) * (s.sx0 - w.xst[n]));
    h.f[AZ - 1] = static_cast<float>(180.0 / PI * std::atan2(0.0f, w.xst[n] - s.sx0));
    h.f[BAZ - 1] = static_cast<float>(180.0 / PI * std::atan2(0.0f, s.sx0 - w.xst[n]));
    return h;
}

void emit(const State& s, SacHeader h, const char* cmp, int idep, const std::vector<float>& d, float cmpinc = -12345.0f,
          float cmpaz = -12345.0f) {
    h.k[KCMPNM] = cmp;
    h.i[IDEP - 71] = idep;
    if (cmpinc != -12345.0f) h.f[CMPINC - 1] = cmpinc;
    if (cmpaz != -12345.0f) h.f[CMPAZ - 1] = cmpaz;
    sac_write(s.odir + "/wav/" + s.title + ".psv." + h.k[KSTNM] + "." + cmp + ".sac", h, d);
}
}  // namespace

void wav_write(const State& s, const Waveforms& w) {
    if (w.nst == 0) return;
    const float mag = moment_magnitude(s.M0);
    for (int n = 0; n < w.nst; ++n) {
        const SacHeader h = base_header(s, w, n, mag);
        if (w.sw_v) {
            emit(s, h, "Vx", 7, w.vel[n][0], 90.0f, 0.0f + s.phi);
            emit(s, h, "Vz", 7, w.vel[n][1], 90.0f, 0.0f);
        }
        if (w.sw_u) {
            emit(s, h, "Ux", 6, w.disp[n][0], 90.0f, 0.0f + s.phi);
            emit(s, h, "Uz", 6, w.disp[n][1], 90.0f, 0.0f);
        }
        if (w.sw_stress) {
            emit(s, h, "Sxx", 5, w.stress[n][0]);
            emit(s, h, "Szz", 5, w.stress[n][1]);
            emit(s, h, "Sxz", 5, w.stress[n][2]);
        }
        if (w.sw_strain) {
            emit(s, h, "Exx", 5, w.strain[n][0]);
            emit(s, h, "Ezz", 5, w.strain[n][1]);
            emit(s, h, "Exz", 5, w.strain[n][2]);
        }
    }
}

}  // namespace swpc
