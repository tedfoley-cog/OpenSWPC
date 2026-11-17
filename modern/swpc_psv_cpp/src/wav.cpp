#include "wav.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "fdtool.hpp"
#include "readini.hpp"
#include "geomap.hpp"
#include "sac.hpp"

namespace swpc {

namespace {
void set_stinfo(State& s, Waveforms& w, const std::string& fn, const std::string& st_format) {
    std::ifstream in(fn);
    if (!in) {
        w.nst = 0;
        return;
    }
    std::string line;
    const float dx = static_cast<float>(s.dx), dz = static_cast<float>(s.dz);
    while (std::getline(in, line)) {
        size_t a = line.find_first_not_of(" \t\r");
        if (a == std::string::npos || line[a] == '#') continue;
        std::string t = line.substr(a);
        for (auto& c : t)
            if (c == ',') c = ' ';
        std::istringstream ss(t);
        std::string a1, a2, a3, stnm, zsw;
        if (!(ss >> a1 >> a2 >> a3 >> stnm >> zsw)) throw std::runtime_error("bad station line: " + line);
        float c1, c2, zst;
        try {
            c1 = static_cast<float>(parse_fortran_real(a1));
            c2 = static_cast<float>(parse_fortran_real(a2));
            zst = static_cast<float>(parse_fortran_real(a3));
        } catch (const std::exception&) {
            throw std::runtime_error("bad station line: " + line);
        }
        for (char q : {'\'', '"'}) {
            if (stnm.size() >= 2 && stnm.front() == q && stnm.back() == q) stnm = stnm.substr(1, stnm.size() - 2);
            if (zsw.size() >= 2 && zsw.front() == q && zsw.back() == q) zsw = zsw.substr(1, zsw.size() - 2);
        }
        if (stnm.size() > 8) stnm = stnm.substr(0, 8);
        if (zsw.size() > 3) zsw = zsw.substr(0, 3);
        float xst, stlo, stla, rdum;
        if (st_format == "xy") {
            xst = c1;
            geomap_c2g(xst, 0.0f, s.clon, s.clat, s.phi, stlo, stla);
        } else if (st_format == "ll") {
            stlo = c1; stla = c2;
            if (!(-360.f <= stlo && stlo <= 360.f && -90.f <= stla && stla <= 90.f))
                throw std::runtime_error("station lon/lat out of range: " + line);
            geomap_g2c(stlo, stla, s.clon, s.clat, s.phi, xst, rdum);
        } else {
            throw std::runtime_error("unknown st_format: " + st_format);
        }
        const int ist = x2i(xst, s.xbeg, dx);
        int kst = z2k(zst, s.zbeg, dz);
        if (!(i2x(1, s.xbeg, dx) < xst && xst < i2x(s.nx, s.xbeg, dx) && s.kbeg < kst && kst < s.kend)) continue;
        if (!(s.ibeg <= ist && ist <= s.iend)) continue;
        if (zsw == "dep") kst = z2k(zst, s.zbeg, dz);
        else if (zsw == "fsb") kst = s.kfs(ist) + 1;
        else if (zsw == "obb") kst = s.kob(ist) + 1;
        else if (zsw == "oba") kst = s.kob(ist) - 1;
        else if (zsw.size() == 3 && zsw[0] == 'b' && zsw[1] == 'd' && zsw[2] >= '0' && zsw[2] <= '9')
            kst = z2k(s.bddep(zsw[2] - '0', ist), s.zbeg, dz);
        else kst = z2k(zst, s.zbeg, dz);
        if (kst > s.kend) kst = s.kend - 1;
        if (kst < s.kbeg) kst = s.kbeg + 1;
        w.xst.push_back(xst); w.zst.push_back(zst); w.stnm.push_back(stnm);
        w.stlo.push_back(stlo); w.stla.push_back(stla); w.ist.push_back(ist); w.kst.push_back(kst);
    }
    w.nst = static_cast<int>(w.ist.size());
}

void alloc3(std::vector<std::vector<std::vector<float>>>& a, int nst, int nc, int ntw) {
    a.assign(nst, std::vector<std::vector<float>>(nc, std::vector<float>(ntw, 0.0f)));
}
}  // namespace

void wav_setup(State& s, Waveforms& w, const IniFile& prm) {
    w.ntdec_w = prm.get_int("ntdec_w", 10);
    w.sw_v = prm.get_bool("sw_wav_v", false);
    w.sw_u = prm.get_bool("sw_wav_u", false);
    w.sw_stress = prm.get_bool("sw_wav_stress", false);
    w.sw_strain = prm.get_bool("sw_wav_strain", false);
    w.wav_format = prm.get_string("wav_format", "sac");
    const std::string st_format = prm.get_string("st_format", "xy");
    const std::string fn_stloc = prm.get_string("fn_stloc", "");
    w.ntdec_w_prg = prm.get_int("ntdec_w_prg", 0);
    if (!(w.sw_v || w.sw_u || w.sw_stress || w.sw_strain)) {
        w.nst = 0;
        return;
    }
    if (w.wav_format != "sac")
        throw std::runtime_error("wav_format '" + w.wav_format + "' is outside the C++ port scope (sac only)");
    w.ntw = static_cast<int>(std::floor(static_cast<float>(s.nt - 1) / static_cast<float>(w.ntdec_w) + 1.0f));
    w.r40x = 9.0 / 8.0 / s.dx;
    w.r40z = 9.0 / 8.0 / s.dz;
    w.r41x = 1.0 / 24.0 / s.dx;
    w.r41z = 1.0 / 24.0 / s.dz;
    set_stinfo(s, w, fn_stloc, st_format);
    std::error_code ec;
    std::filesystem::create_directories(s.odir + "/wav", ec);
    if (w.sw_v) alloc3(w.vel, w.nst, 2, w.ntw);
    if (w.sw_u) alloc3(w.disp, w.nst, 2, w.ntw);
    if (w.sw_stress) alloc3(w.stress, w.nst, 3, w.ntw);
    if (w.sw_strain) alloc3(w.strain, w.nst, 3, w.ntw);
    w.ux.assign(w.nst, 0.0f); w.uz.assign(w.nst, 0.0f);
    w.exx.assign(w.nst, 0.0f); w.ezz.assign(w.nst, 0.0f); w.exz.assign(w.nst, 0.0f);
}

void wav_store(State& s, Waveforms& w, int it) {
    if (w.nst == 0) return;
    const float dt = s.dt;
    if (w.sw_u) {
        for (int n = 0; n < w.nst; ++n) {
            const int i = w.ist[n], k = w.kst[n];
            w.ux[n] = w.ux[n] + (s.Vx(k, i) + s.Vx(k, i - 1)) * 0.5f * dt;
            w.uz[n] = w.uz[n] - (s.Vz(k, i) + s.Vz(k - 1, i)) * 0.5f * dt;
        }
    }
    if (w.sw_strain) {
        const real_mp r40x = w.r40x, r40z = w.r40z, r41x = w.r41x, r41z = w.r41z;
        for (int n = 0; n < w.nst; ++n) {
            const int i = w.ist[n], k = w.kst[n];
            auto& Vx = s.Vx;
            auto& Vz = s.Vz;
            const real_mp dxVx = (Vx(k, i) - Vx(k, i - 1)) * r40x - (Vx(k, i + 1) - Vx(k, i - 2)) * r41x;
            const real_mp dzVz = (Vz(k, i) - Vz(k - 1, i)) * r40z - (Vz(k + 1, i) - Vz(k - 2, i)) * r41z;
            const real_mp dxVz = ((Vz(k, i + 1) - Vz(k, i)) * r40x - (Vz(k, i + 2) - Vz(k, i - 1)) * r41x +
                                  (Vz(k - 1, i + 1) - Vz(k - 1, i)) * r40x - (Vz(k - 1, i + 2) - Vz(k - 1, i - 1)) * r41x +
                                  (Vz(k, i) - Vz(k, i - 1)) * r40x - (Vz(k, i + 1) - Vz(k, i - 2)) * r41x +
                                  (Vz(k - 1, i) - Vz(k - 1, i - 1)) * r40x - (Vz(k - 1, i + 1) - Vz(k - 1, i - 2)) * r41x) /
                                 4.0;
            const real_mp dzVx = ((Vx(k + 1, i) - Vx(k, i)) * r40z - (Vx(k + 2, i) - Vx(k - 1, i)) * r41z +
                                  (Vx(k + 1, i - 1) - Vx(k, i - 1)) * r40z - (Vx(k + 2, i - 1) - Vx(k - 1, i - 1)) * r41z +
                                  (Vx(k, i) - Vx(k - 1, i)) * r40z - (Vx(k + 1, i) - Vx(k - 2, i)) * r41z +
                                  (Vx(k, i - 1) - Vx(k - 1, i - 1)) * r40z - (Vx(k + 1, i - 1) - Vx(k - 2, i - 1)) * r41z) /
                                 4.0;
            w.exx[n] = w.exx[n] + static_cast<float>(dxVx) * dt;
            w.ezz[n] = w.ezz[n] + static_cast<float>(dzVz) * dt;
            w.exz[n] = w.exz[n] + static_cast<float>(dxVz + dzVx) / 2.0f * dt;
        }
    }
    if ((it - 1) % w.ntdec_w == 0) {
        const int itw = (it - 1) / w.ntdec_w;
        const float sc = s.M0 * s.UC;
        for (int n = 0; n < w.nst; ++n) {
            const int i = w.ist[n], k = w.kst[n];
            if (w.sw_v) {
                w.vel[n][0][itw] = (s.Vx(k, i) + s.Vx(k, i - 1)) / 2.0f * s.M0 * s.UC * 1e9f;
                w.vel[n][1][itw] = -(s.Vz(k, i) + s.Vz(k - 1, i)) / 2.0f * s.M0 * s.UC * 1e9f;
            }
            if (w.sw_u) {
                w.disp[n][0][itw] = w.ux[n] * s.M0 * s.UC * 1e9f;
                w.disp[n][1][itw] = w.uz[n] * s.M0 * s.UC * 1e9f;
            }
            if (w.sw_stress) {
                w.stress[n][0][itw] = static_cast<float>(s.Sxx(k, i)) * s.M0 * s.UC * 1e6f;
                w.stress[n][1][itw] = static_cast<float>(s.Szz(k, i)) * s.M0 * s.UC * 1e6f;
                w.stress[n][2][itw] = static_cast<float>(s.Sxz(k, i) + s.Sxz(k, i - 1) + s.Sxz(k - 1, i) +
                                                         s.Sxz(k - 1, i - 1)) / 4.0f * s.M0 * s.UC * 1e6f;
            }
            if (w.sw_strain) {
                w.strain[n][0][itw] = w.exx[n] * s.M0 * s.UC * 1e-3f;
                w.strain[n][1][itw] = w.ezz[n] * s.M0 * s.UC * 1e-3f;
                w.strain[n][2][itw] = w.exz[n] * s.M0 * s.UC * 1e-3f;
            }
        }
        (void)sc;
    }
    if (w.ntdec_w_prg > 0 && (it - 1) % w.ntdec_w_prg == 0) wav_write(s, w);
}

// wav_write is implemented in wav_write.cpp (SAC header population).

}  // namespace swpc
