// swpc_psv (C++17): port of src/swpc_psv/main.f90. Same input (-i file.inf) and SAC output layout.
#include <chrono>
#include <cstdio>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <string>

#include "dump.hpp"
#include "fdtool.hpp"
#include "global.hpp"
#include "kernel.hpp"
#include "medium.hpp"
#include "readini.hpp"
#include "source.hpp"
#include "state.hpp"
#include "wav.hpp"

using namespace swpc;

namespace {
void usage() { std::fprintf(stderr, "usage: swpc_psv_cpp -i <parameter file> [-q]\n"); }

void report_setup(const State& s, const Source& src, const Waveforms& w, bool quiet) {
    float c = 0, r = 0;
    fdm_cond_stability(static_cast<float>(s.dx), 1e10f, static_cast<float>(s.dz), s.vmax, s.dt, c);
    fdm_cond_wavelength(static_cast<float>(s.dx), -1.0f, static_cast<float>(s.dz), s.vmin, s.fmax, r);
    if (c > 1.0f) throw std::runtime_error("CFL stability condition violated (c = " + std::to_string(c) + " > 1)");
    if (quiet) return;
    std::printf("swpc_psv_cpp: title=%s nx=%d nz=%d nt=%d dx=%g dz=%g dt=%g\n", s.title.c_str(), s.nx, s.nz, s.nt,
                s.dx, s.dz, s.dt);
    std::printf("  vmin=%g vmax=%g fmax=%g stability=%g wavelength_ppw=%g\n", s.vmin, s.vmax, s.fmax, c, r);
    std::printf("  abc_type=%s na=%d nsrc=%d nst=%d ntw=%d\n", s.abc_type.c_str(), s.na, src.nsrc, w.nst, w.ntw);
}
}  // namespace

int main(int argc, char** argv) {
    std::string fn_prm;
    bool quiet = false;
    for (int a = 1; a < argc; ++a) {
        if (!std::strcmp(argv[a], "-i") && a + 1 < argc) fn_prm = argv[++a];
        else if (!std::strcmp(argv[a], "-q")) quiet = true;
        else { usage(); return 2; }
    }
    if (fn_prm.empty()) { usage(); return 2; }
    try {
        const IniFile prm(fn_prm);
        State s;
        Source src;
        Waveforms w;
        Absorber ab;
        global_setup(s, prm);
        global_setup2(s);
        medium_setup(s, prm);
        kernel_setup(s);
        source_setup(s, src, prm);
        absorb_setup(s, ab);
        dump_medium(s);
        wav_setup(s, w, prm);
        report_setup(s, src, w, quiet);

        const auto t0 = std::chrono::steady_clock::now();
        for (int it = 1; it <= s.nt; ++it) {
            wav_store(s, w, it);
            kernel_update_stress(s);
            absorb_update_stress(s, ab);
            source_stressglut(s, src, it);
            kernel_update_vel(s);
            source_bodyforce(s, src, it);
            absorb_update_vel(s, ab);
            dump_state(s, it);
        }
        wav_write(s, w);
        if (!quiet) {
            const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::printf("  time loop: %.3f s (%.3g cell-updates/s)\n", sec,
                        static_cast<double>(s.nx) * s.nz * s.nt / (sec > 0 ? sec : 1));
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "ERROR: %s\n", e.what());
        return 1;
    }
    return 0;
}
