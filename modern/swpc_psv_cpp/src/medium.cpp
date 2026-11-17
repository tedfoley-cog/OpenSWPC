#include "medium.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "fdtool.hpp"

namespace swpc {

namespace {
constexpr double R_EARTH = 6371.0;
bool g_use_munk = false;

std::string strip(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

void flattening(const State& s, bool earth_flattening, Array1D<float>& zs, Array1D<float>& Cv) {
    zs.allocate(s.kbeg_m, s.kend_m);
    Cv.allocate(s.kbeg_m, s.kend_m);
    for (int k = s.kbeg_m; k <= s.kend_m; ++k) {
        if (earth_flattening) {
            zs(k) = static_cast<float>(R_EARTH - R_EARTH * std::exp(-s.zc(k) / R_EARTH));
            Cv(k) = static_cast<float>(std::exp(s.zc(k) / R_EARTH));
        } else {
            zs(k) = s.zc(k);
            Cv(k) = 1.0f;
        }
    }
}

void set_cell(State& s, int k, int i, float rho1, float vp1, float vs1, float qp1, float qs1) {
    s.rho(k, i) = rho1;
    s.mu(k, i) = rho1 * vs1 * vs1;
    s.lam(k, i) = rho1 * (vp1 * vp1 - 2 * vs1 * vs1);
    s.taup(k, i) = qp1;  // holds Q until converted to tau below
    s.taus(k, i) = qs1;
}

void vmodel_uni(State& s, const IniFile& prm) {
    const float vp0 = prm.get_float("vp0", 5.0f);
    const float vs0 = prm.get_float("vs0", vp0 / std::sqrt(3.0f));
    const float rho0 = prm.get_float("rho0", 2.7f);
    const float qp0 = prm.get_float("qp0", 1000000.0f);
    const float qs0 = prm.get_float("qs0", 1000000.0f);
    const float topo0 = prm.get_float("topo0", 0.0f);
    g_use_munk = prm.get_bool("munk_profile", false);
    const bool ef = prm.get_bool("earth_flattening", false);
    Array1D<float> zs, Cv;
    flattening(s, ef, zs, Cv);
    for (int i = s.ibeg_m; i <= s.iend_m; ++i) {
        if (!s.fullspace_mode) s.bddep(0, i) = topo0;
        for (int k = s.kbeg_m; k <= s.kend_m; ++k) {
            if (s.fullspace_mode || zs(k) > topo0) {
                set_cell(s, k, i, rho0, Cv(k) * vp0, Cv(k) * vs0, qp0, qs0);
            } else if (s.zc(k) > 0.0f) {
                set_cell(s, k, i, 1.0f, Cv(k) * seawater_vel(zs(k), g_use_munk), 0.0f, 1000000.0f, 1000000.0f);
            } else {
                set_cell(s, k, i, 0.001f, 0.0f, 0.0f, 10.0f, 10.0f);
            }
        }
        for (int j = 1; j <= NBD; ++j) s.bddep(j, i) = -9999;
    }
}

void vmodel_lhm(State& s, const IniFile& prm) {
    const std::string fn = prm.get_string("fn_lhm", "");
    const float vcut = prm.get_float("vcut", 0.0f);
    std::ifstream in(fn);
    if (!in) throw std::runtime_error("vmodel_lhm: file not found: " + fn);
    g_use_munk = prm.get_bool("munk_profile", false);
    const bool ef = prm.get_bool("earth_flattening", false);
    Array1D<float> zs, Cv;
    flattening(s, ef, zs, Cv);

    std::vector<float> depth, rho0, vp0, vs0, qp0, qs0;
    std::string line;
    while (std::getline(in, line)) {
        std::string a = strip(line);
        if (a.empty() || a[0] == '#') continue;
        for (auto& c : a)
            if (c == ',') c = ' ';
        std::istringstream ss(a);
        float d, r, p, v, qp, qs;
        if (!(ss >> d >> r >> p >> v >> qp >> qs)) throw std::runtime_error("vmodel_lhm: bad line: " + line);
        depth.push_back(d); rho0.push_back(r); vp0.push_back(p); vs0.push_back(v); qp0.push_back(qp); qs0.push_back(qs);
    }
    const int nlayer = static_cast<int>(depth.size());
    if (nlayer == 0) throw std::runtime_error("vmodel_lhm: no layers in " + fn);
    for (int l = nlayer - 2; l >= 0; --l) {
        if ((vp0[l] < vcut || vs0[l] < vcut) && (vp0[l] > 0 && vs0[l] > 0)) {
            vp0[l] = vp0[l + 1]; vs0[l] = vs0[l + 1]; rho0[l] = rho0[l + 1];
            qp0[l] = qp0[l + 1]; qs0[l] = qs0[l + 1];
        }
    }
    for (int i = s.ibeg_m; i <= s.iend_m; ++i) s.bddep(0, i) = depth[0];
    float rho1 = 0.001f, vp1 = 0.0f, vs1 = 0.0f, qp1 = 10.0f, qs1 = 10.0f;
    for (int k = s.kbeg_m; k <= s.kend_m; ++k) {
        if (zs(k) < depth[0]) {
            if (zs(k) < 0.0f) {
                rho1 = 0.001f; vp1 = 0.0f; vs1 = 0.0f; qp1 = 10.0f; qs1 = 10.0f;
            } else {
                rho1 = 1.0f; vp1 = Cv(k) * seawater_vel(zs(k), g_use_munk); vs1 = 0.0f;
                qp1 = 1000000.0f; qs1 = 1000000.0f;
            }
        } else {
            for (int l = 0; l < nlayer; ++l) {
                if (zs(k) >= depth[l]) {
                    rho1 = rho0[l]; vp1 = Cv(k) * vp0[l]; vs1 = Cv(k) * vs0[l]; qp1 = qp0[l]; qs1 = qs0[l];
                }
            }
        }
        for (int i = s.ibeg_m; i <= s.iend_m; ++i) set_cell(s, k, i, rho1, vp1, vs1, qp1, qs1);
    }
    for (int i = s.ibeg_m; i <= s.iend_m; ++i)
        for (int j = 1; j <= NBD; ++j) s.bddep(j, i) = -9999;
}

template <class F>
void for_all_media(State& s, F&& f) {
    f(s.rho); f(s.lam); f(s.mu); f(s.taup); f(s.taus);
}

void stabilize_absorber(State& s) {
    const float V_DYNAMIC_RANGE = 0.4f;
    const int LV_THICK = 20;
    const float vmin_pml = s.vmax * V_DYNAMIC_RANGE;
    for (const Box& b : s.box) {
        if (b.ncell == 0) continue;
        for (int i = b.ib; i <= b.ie; ++i) {
            int k = b.kb + 1;
            while (k <= b.ke) {
                if (s.lam(k, i) < s.lam(k - 1, i) || s.mu(k, i) < s.mu(k - 1, i)) {
                    int k2 = k + 1;
                    for (; k2 <= s.kend; ++k2)
                        if (s.lam(k2, i) > s.lam(k2 - 1, i) || s.mu(k2, i) > s.mu(k2 - 1, i)) break;
                    if (k2 - k <= LV_THICK) {
                        for_all_media(s, [&](Array2D<float>& a) {
                            for (int kk = k; kk <= k2 - 1; ++kk) a(kk, i) = a(k - 1, i);
                        });
                        k = k2 - 1;
                    }
                }
                ++k;
            }
        }
    }
    for (const Box& b : s.box) {
        if (b.ncell == 0) continue;
        for (int i = b.ib; i <= b.ie; ++i) {
            for (int k = b.kb; k <= b.ke; ++k) {
                float vp = std::sqrt((s.lam(k, i) + 2 * s.mu(k, i)) / s.rho(k, i));
                float vs = std::sqrt(s.mu(k, i) / s.rho(k, i));
                (void)vp;
                if (vs < FLT_EPSILON) continue;
                const float gamma = std::sqrt(3.0f);
                if (vs < vmin_pml) {
                    vs = vmin_pml;
                    vp = vs * gamma;
                    s.lam(k, i) = s.rho(k, i) * (vp * vp - 2 * vs * vs);
                    s.mu(k, i) = s.rho(k, i) * (vs * vs);
                }
            }
        }
    }
}
}  // namespace

float seawater_vel(float z, bool use_munk) {
    const double epsil = use_munk ? 0.00737 : 0.0;
    const double zc = 1300.0;
    const double zb = 2 * (static_cast<double>(z) * 1000. - zc) / zc;
    return static_cast<float>(1.5 * (1.0 + epsil * (zb - 1.0 + std::exp(-zb))));
}

void medium_setup(State& s, const IniFile& prm) {
    const int kb = s.kbeg_m, ke = s.kend_m, ib = s.ibeg_m, ie = s.iend_m;
    for_all_media(s, [&](Array2D<float>& a) { a.allocate(kb, ke, ib, ie, 0.0f); });
    for (auto* a : {&s.kfs, &s.kob, &s.kfs_top, &s.kfs_bot, &s.kob_top, &s.kob_bot}) a->allocate(ib, ie, 0);
    s.bddep.allocate(0, NBD, ib, ie, 0.0f);
    s.ts.assign(NM, 0.0f);

    float fq_min, fq_max, fq_ref;
    if (s.benchmark_mode) {
        fq_min = 0.05f; fq_max = 5.0f; fq_ref = 1.0f;
        for (int k = kb; k <= ke; ++k)
            for (int i = ib; i <= ie; ++i) {
                if (s.zc(k) < 0.0f) {
                    s.rho(k, i) = 0.001f; s.mu(k, i) = 0.0f; s.lam(k, i) = 0.0f;
                } else {
                    s.rho(k, i) = 2.7f; s.mu(k, i) = 2.7f * 3.5f * 3.5f; s.lam(k, i) = 2.7f * 3.5f * 3.5f;
                }
                s.taup(k, i) = 1e10f; s.taus(k, i) = 1e10f;
            }
    } else {
        fq_min = prm.get_float("fq_min", 0.05f);
        fq_max = prm.get_float("fq_max", 5.00f);
        fq_ref = prm.get_float("fq_ref", 1.00f);
        const std::string vmodel_type = prm.get_string("vmodel_type", "uni");
        if (vmodel_type == "uni") vmodel_uni(s, prm);
        else if (vmodel_type == "lhm") vmodel_lhm(s, prm);
        else throw std::runtime_error("vmodel_type '" + vmodel_type + "' is not supported by the C++ port (supported: uni, lhm)");
    }

    const int na = s.na, nx = s.nx, nz = s.nz;
    for_all_media(s, [&](Array2D<float>& a) {
        for (int i = ib; i <= na; ++i)
            for (int k = kb; k <= ke; ++k) a(k, i) = a(k, na + 1);
        for (int i = nx - na + 1; i <= ie; ++i)
            for (int k = kb; k <= ke; ++k) a(k, i) = a(k, nx - na);
        for (int i = ib; i <= ie; ++i)
            for (int k = nz - na + 1; k <= ke; ++k) a(k, i) = a(nz - na, i);
    });

    visco_set_relaxtime(NM, s.ts, fq_min, fq_max);
    float zeta = visco_constq_zeta(NM, fq_min, fq_max, s.ts);
    if (s.benchmark_mode) zeta = 0.0f;
    for (int i = ib; i <= ie; ++i)
        for (int k = kb; k <= ke; ++k) {
            s.taup(k, i) = NM * zeta / s.taup(k, i);
            s.taus(k, i) = NM * zeta / s.taus(k, i);
        }

    // relaxed medium
    for (int i = ib; i <= ie; ++i)
        for (int k = kb; k <= ke; ++k) {
            const float rho_beta2 = s.mu(k, i);
            const float rho_alpha2 = s.lam(k, i) + 2 * s.mu(k, i);
            const float cs = visco_chi(NM, s.ts, s.taus(k, i), fq_ref);
            const float cp = visco_chi(NM, s.ts, s.taup(k, i), fq_ref);
            s.mu(k, i) = rho_beta2 / (cs * cs);
            s.lam(k, i) = rho_alpha2 / (cp * cp) - 2 * s.mu(k, i);
        }

    // surface detection
    const float epsl = FLT_EPSILON;
    s.kfs.fill(s.kbeg - 1);
    s.kob.fill(s.kbeg - 1);
    for (int i = s.ibeg - 1; i <= s.iend + 2; ++i)
        for (int k = s.kbeg; k <= s.kend - 1; ++k) {
            if (std::fabs(s.mu(k, i)) < epsl && std::fabs(s.mu(k + 1, i)) > epsl) s.kob(i) = k;
            if (std::fabs(s.lam(k, i)) < epsl && std::fabs(s.lam(k + 1, i)) > epsl) s.kfs(i) = k;
        }
    for (int i = s.ibeg; i <= s.iend; ++i) {
        int fmin = s.kfs(i - 2), fmx = s.kfs(i - 2), omin = s.kob(i - 2), omx = s.kob(i - 2);
        for (int ii = i - 2; ii <= i + 3; ++ii) {
            fmin = std::min(fmin, s.kfs(ii)); fmx = std::max(fmx, s.kfs(ii));
            omin = std::min(omin, s.kob(ii)); omx = std::max(omx, s.kob(ii));
        }
        s.kfs_top(i) = std::max(fmin - 2, s.kbeg);
        s.kfs_bot(i) = std::min(fmx + 2, s.kend);
        s.kob_top(i) = std::max(omin - 2, s.kbeg);
        s.kob_bot(i) = std::min(omx + 2, s.kend);
    }
    if (s.fullspace_mode) {
        for (int i = ib; i <= ie; ++i) {
            s.kfs_bot(i) = s.kfs_top(i) - 1;
            s.kob_bot(i) = s.kob_top(i) - 1;
        }
    }

    // velocity min/max
    float vmax1 = -1, vmin1 = 1e30f;
    for (int i = s.ibeg; i <= s.iend; ++i)
        for (int k = s.kfs(i) + 1; k <= s.kend; ++k) {
            const float vp = std::sqrt((s.lam(k, i) + 2 * s.mu(k, i)) / s.rho(k, i));
            vmax1 = std::max(vmax1, vp);
            const float vs = std::sqrt(s.mu(k, i) / s.rho(k, i));
            if (vs < FLT_EPSILON) continue;
            vmin1 = std::min(vmin1, vs);
        }
    s.vmax = vmax1;
    s.vmin = vmin1;

    if (prm.get_bool("stabilize_pml", false)) stabilize_absorber(s);
}

}  // namespace swpc
