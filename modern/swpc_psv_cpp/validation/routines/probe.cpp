// C++ counterpart of probe.f90: identical inputs, identical bit-pattern output.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "fdtool.hpp"
#include "geomap.hpp"
#include "kernel.hpp"
#include "medium.hpp"
#include "state.hpp"

using namespace swpc;

static unsigned hx(float v) {
    std::uint32_t u;
    std::memcpy(&u, &v, 4);
    return u;
}

int main() {
    const char* stf[7] = {"kupper", "texp", "cosine", "boxcar", "triangle", "herrmann", "asymcos"};
    for (int j = 0; j < 7; ++j)
        for (int i = 0; i <= 400; ++i) {
            const float t = -0.2f + i * 0.005f;
            const float prm[3] = {0.1f, 0.7f, 1.1f};
            std::printf("stf %s%5d %08X\n", stf[j], i, hx(momentrate(t, stf[j], prm)));
        }
    for (int i = 0; i <= 36; ++i)
        for (int j = 0; j <= 9; ++j)
            for (int k = 0; k <= 12; ++k) {
                float m[6];
                sdr2moment(i * 10.0f - 180.0f, j * 10.0f, k * 30.0f - 180.0f, m[0], m[1], m[2], m[3], m[4], m[5]);
                std::printf("sdr%4d%4d%4d", i, j, k);
                for (float v : m) std::printf(" %08X", hx(v));
                std::printf("\n");
            }
    for (int i = 0; i <= 80; ++i)
        std::printf("mw%4d %08X %08X\n", i, hx(seismic_moment(i * 0.1f)),
                    hx(moment_magnitude(std::pow(10.0f, 10 + i * 0.137f))));
    for (int i = -50; i <= 50; ++i) {
        const float x = i * 0.0731f;
        std::printf("grid%4d%8d%8d %08X %08X\n", i, x2i(x, -1.3f, 0.025f), z2k(x, -0.7f, 0.0125f),
                    hx(i2x(i, -1.3f, 0.025f)), hx(k2z(i, -0.7f, 0.0125f)));
    }
    for (int i = 1; i <= 6; ++i) {
        std::vector<float> ts(3);
        visco_set_relaxtime(3, ts, 0.01f * i, 2.0f * i);
        const float zeta = visco_constq_zeta(3, 0.01f * i, 2.0f * i, ts);
        std::printf("visco%4d %08X %08X %08X %08X %08X\n", i, hx(ts[0]), hx(ts[1]), hx(ts[2]), hx(zeta),
                    hx(visco_chi(3, ts, zeta, 1.0f)));
    }
    for (int i = 0; i <= 40; ++i) {
        float lon = 139.0f + i * 0.05f, lat = 35.0f + i * 0.03f, x, y;
        geomap_g2c(lon, lat, 139.7604f, 35.7182f, 12.5f, x, y);
        geomap_c2g(x, y, 139.7604f, 35.7182f, 12.5f, lon, lat);
        std::printf("geo%4d %08X %08X %08X %08X\n", i, hx(x), hx(y), hx(lon), hx(lat));
    }
    for (int i = 0; i <= 60; ++i) std::printf("munk%4d %08X\n", i, hx(seawater_vel(i * 0.1f, true)));
    for (int j = 1; j <= 3; ++j) {
        State s;
        s.na = 8 * j; s.fcut = 0.5f * j; s.dt = 0.004f * j;
        for (int i = -5; i <= 125; ++i) {
            float g[4];
            damping_profile(s, i * 0.1f, s.na * 0.1f, 0.0f, 12.0f, g);
            std::printf("pml%4d%4d %08X %08X %08X %08X\n", j, i, hx(g[0]), hx(g[1]), hx(g[2]), hx(g[3]));
        }
    }
    return 0;
}
