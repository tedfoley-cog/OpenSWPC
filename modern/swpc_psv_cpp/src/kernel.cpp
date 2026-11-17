#include "kernel.hpp"

#include <cfloat>
#include <cmath>
#include <stdexcept>

#include "fdtool.hpp"

namespace swpc {

namespace {
inline int isign_of(int v) { return v >= 0 ? 1 : -1; }  // Fortran sign(1, v)

inline int order_sign(const State& s, int k, int i) {
    return isign_of(std::max((k - s.kfs_top(i)) * (s.kfs_bot(i) - k), (k - s.kob_top(i)) * (s.kob_bot(i) - k)));
}

inline float harmonic_mu(float nnn, float pnn, float npn, float ppn) {
    const float epsl = FLT_EPSILON;
    return 4 * nnn * pnn * npn * ppn /
           (nnn * pnn * npn + nnn * pnn * ppn + nnn * npn * ppn + pnn * npn * ppn + epsl);
}
}  // namespace

FdCoef make_fd_coef(const State& s) {
    FdCoef c;
    c.rc40x = 17.0 / 16.0 / s.dx;
    c.rc40z = 17.0 / 16.0 / s.dz;
    c.rc41x = 1.0 / 48.0 / s.dx;
    c.rc41z = 1.0 / 48.0 / s.dz;
    c.rd40x = -1.0 / 16.0 / s.dx;
    c.rd40z = -1.0 / 16.0 / s.dz;
    c.rd41x = -1.0 / 48.0 / s.dx;
    c.rd41z = -1.0 / 48.0 / s.dz;
    return c;
}

void kernel_setup(State& s) {
    const int kb = s.kbeg_m, ke = s.kend_m, ib = s.ibeg_m, ie = s.iend_m;
    s.Vx.allocate(kb, ke, ib, ie, 0.0f);
    s.Vz.allocate(kb, ke, ib, ie, 0.0f);
    s.Sxx.allocate(kb, ke, ib, ie, 0.0);
    s.Szz.allocate(kb, ke, ib, ie, 0.0);
    s.Sxz.allocate(kb, ke, ib, ie, 0.0);
    s.Rxx.allocate(NM, s.kbeg, s.kend, s.ibeg, s.iend, 0.0f);
    s.Rzz.allocate(NM, s.kbeg, s.kend, s.ibeg, s.iend, 0.0f);
    s.Rxz.allocate(NM, s.kbeg, s.kend, s.ibeg, s.iend, 0.0f);
    s.c1.assign(NM, 0.0f);
    s.c2.assign(NM, 0.0f);
    s.d1.assign(NM, 0.0f);
    const float dt = s.dt;
    float d2sum = 0.0f;
    for (int m = 0; m < NM; ++m) {
        const float t = s.ts[m];
        s.c1[m] = (2 * t - dt) / (2 * t + dt);
        s.c2[m] = 2.0f / (2 * t + dt) / NM;
        d2sum += dt / (2 * t - dt);
        s.d1[m] = 2 * t / (2 * t - dt);
    }
    s.d2 = d2sum / NM;
}

void kernel_update_vel(State& s) {
    const FdCoef c = make_fd_coef(s);
    const float dt = s.dt;
#pragma omp parallel for schedule(static, 1)
    for (int i = s.ibeg_k; i <= s.iend_k; ++i) {
        for (int k = s.kbeg_k; k <= s.kend_k; ++k) {
            const int isign = order_sign(s, k, i);
            const real_mp re40x = c.rc40x + isign * c.rd40x, re41x = c.rc41x + isign * c.rd41x;
            const real_mp re40z = c.rc40z + isign * c.rd40z, re41z = c.rc41z + isign * c.rd41z;
            const real_mp dxSxx = (s.Sxx(k, i + 1) - s.Sxx(k, i)) * re40x - (s.Sxx(k, i + 2) - s.Sxx(k, i - 1)) * re41x;
            const real_mp dzSzz = (s.Szz(k + 1, i) - s.Szz(k, i)) * re40z - (s.Szz(k + 2, i) - s.Szz(k - 1, i)) * re41z;
            const real_mp dxSxz = (s.Sxz(k, i) - s.Sxz(k, i - 1)) * re40x - (s.Sxz(k, i + 1) - s.Sxz(k, i - 2)) * re41x;
            const real_mp dzSxz = (s.Sxz(k, i) - s.Sxz(k - 1, i)) * re40z - (s.Sxz(k + 1, i) - s.Sxz(k - 2, i)) * re41z;
            const float bx = 2.0f / (s.rho(k, i) + s.rho(k, i + 1));
            const float bz = 2.0f / (s.rho(k, i) + s.rho(k + 1, i));
            s.Vx(k, i) = static_cast<float>(s.Vx(k, i) + bx * (dxSxx + dzSxz) * dt);
            s.Vz(k, i) = static_cast<float>(s.Vz(k, i) + bz * (dxSxz + dzSzz) * dt);
        }
    }
}

void kernel_update_stress(State& s) {
    const FdCoef c = make_fd_coef(s);
    const float dt = s.dt;
    const float d2 = s.d2;
#pragma omp parallel for schedule(static, 1)
    for (int i = s.ibeg_k; i <= s.iend_k; ++i) {
        for (int k = s.kbeg_k; k <= s.kend_k; ++k) {
            const int isign = order_sign(s, k, i);
            const real_mp re40x = c.rc40x + isign * c.rd40x, re41x = c.rc41x + isign * c.rd41x;
            const real_mp re40z = c.rc40z + isign * c.rd40z, re41z = c.rc41z + isign * c.rd41z;
            const real_mp dxVx = (s.Vx(k, i) - s.Vx(k, i - 1)) * re40x - (s.Vx(k, i + 1) - s.Vx(k, i - 2)) * re41x;
            const real_mp dzVz = (s.Vz(k, i) - s.Vz(k - 1, i)) * re40z - (s.Vz(k + 1, i) - s.Vz(k - 2, i)) * re41z;
            const float mu2 = 2 * s.mu(k, i);
            const float lam2mu = s.lam(k, i) + mu2;
            const float taup1 = s.taup(k, i);
            const float taus1 = s.taus(k, i);
            const float d2v2 = static_cast<float>(dxVx + dzVz);
            const float f_Rxx = static_cast<float>(lam2mu * taup1 * d2v2 - mu2 * taus1 * dzVz);
            const float f_Rzz = static_cast<float>(lam2mu * taup1 * d2v2 - mu2 * taus1 * dxVx);
            float Rxx_n = 0.0f, Rzz_n = 0.0f;
            for (int m = 1; m <= NM; ++m) {
                s.Rxx(m, k, i) = s.c1[m - 1] * s.Rxx(m, k, i) - s.c2[m - 1] * f_Rxx * dt;
                s.Rzz(m, k, i) = s.c1[m - 1] * s.Rzz(m, k, i) - s.c2[m - 1] * f_Rzz * dt;
                Rxx_n += s.d1[m - 1] * s.Rxx(m, k, i);
                Rzz_n += s.d1[m - 1] * s.Rzz(m, k, i);
            }
            const float taup_plus1 = 1 + taup1 * (1 + d2);
            const float taus_plus1 = 1 + taus1 * (1 + d2);
            s.Sxx(k, i) += (lam2mu * taup_plus1 * d2v2 - mu2 * taus_plus1 * dzVz + Rxx_n) * dt;
            s.Szz(k, i) += (lam2mu * taup_plus1 * d2v2 - mu2 * taus_plus1 * dxVx + Rzz_n) * dt;
        }
    }
#pragma omp parallel for schedule(static, 1)
    for (int i = s.ibeg_k; i <= s.iend_k; ++i) {
        for (int k = s.kbeg_k; k <= s.kend_k; ++k) {
            const int isign = order_sign(s, k, i);
            const real_mp re40x = c.rc40x + isign * c.rd40x, re41x = c.rc41x + isign * c.rd41x;
            const real_mp re40z = c.rc40z + isign * c.rd40z, re41z = c.rc41z + isign * c.rd41z;
            const real_mp dxVz = (s.Vz(k, i + 1) - s.Vz(k, i)) * re40x - (s.Vz(k, i + 2) - s.Vz(k, i - 1)) * re41x;
            const real_mp dzVx = (s.Vx(k + 1, i) - s.Vx(k, i)) * re40z - (s.Vx(k + 2, i) - s.Vx(k - 1, i)) * re41z;
            const float taus1 = s.taus(k, i);
            const float mu_xz = harmonic_mu(s.mu(k, i), s.mu(k + 1, i), s.mu(k, i + 1), s.mu(k + 1, i + 1));
            const float dxVz_dzVx = static_cast<float>(dxVz + dzVx);
            const float f_Rxz = mu_xz * taus1 * dxVz_dzVx;
            float Rxz_n = 0.0f;
            for (int m = 1; m <= NM; ++m) {
                s.Rxz(m, k, i) = s.c1[m - 1] * s.Rxz(m, k, i) - s.c2[m - 1] * f_Rxz * dt;
                Rxz_n += s.d1[m - 1] * s.Rxz(m, k, i);
            }
            const float taus_plus1 = 1 + taus1 * (1 + d2);
            s.Sxz(k, i) += (mu_xz * taus_plus1 * dxVz_dzVx + Rxz_n) * dt;
        }
    }
}

// ---------------------------------------------------------------- absorbers

void damping_profile(const State& s, float x, float H, float xb, float xe, float g[4]) {
    const int pd = 2;
    const float cp = 6.0f;
    const float R0 = std::pow(10.0f, -(std::log10(static_cast<float>(s.na)) - 1) / std::log10(2.0f) - 3.0f);
    const float d0 = -(1.0f / (2.0f * H)) * (pd + 1) * cp * std::log(R0);
    const float b0 = 7.0f;
    const float a0 = static_cast<float>(PI * s.fcut);
    float xx;
    if (x <= xb + H) xx = (xb + H) - x;
    else if (x >= xe - H) xx = x - (xe - H);
    else xx = 0.0f;
    const float r = std::fabs(xx / H);
    const float d = d0 * (r * r);  // Fortran: d0 * abs(xx/H)**pd, the power is evaluated first
    const float a = a0 * (1.0f - r);
    const float b = 1.0f + (b0 - 1.0f) * (r * r * r);
    const float h = s.dt / 2.0f;
    const float den = 1.0f + h * (a + d / b);
    g[0] = ((1.0f + h * a) / b) / den;
    g[1] = (-1.0f / b) / den;
    g[2] = (1.0f - h * (a + d / b)) / den;
    g[3] = (d / b) / den;
}

namespace {

void pml_setup(State& s, Absorber& a) {
    a.gxc.allocate(1, 4, s.ibeg, s.iend);
    a.gxe.allocate(1, 4, s.ibeg, s.iend);
    a.gzc.allocate(1, 4, s.kbeg, s.kend);
    a.gze.allocate(1, 4, s.kbeg, s.kend);
    const float hx = static_cast<float>(s.na * s.dx);
    const float hz = static_cast<float>(s.na * s.dz);
    float g[4];
    for (int i = s.ibeg; i <= s.iend; ++i) {
        damping_profile(s, s.xc(i), hx, s.xbeg, s.xend, g);
        for (int j = 0; j < 4; ++j) a.gxc(j + 1, i) = g[j];
        damping_profile(s, s.xc(i) + static_cast<float>(s.dx) / 2.0f, hx, s.xbeg, s.xend, g);
        for (int j = 0; j < 4; ++j) a.gxe(j + 1, i) = g[j];
    }
    for (int k = s.kbeg; k <= s.kend; ++k) {
        damping_profile(s, s.zc(k), hz, s.zbeg, s.zend, g);
        for (int j = 0; j < 4; ++j) a.gzc(j + 1, k) = g[j];
        damping_profile(s, s.zc(k) + static_cast<float>(s.dz) / 2.0f, hz, s.zbeg, s.zend, g);
        for (int j = 0; j < 4; ++j) a.gze(j + 1, k) = g[j];
    }
    for (auto* v : {&a.axVx, &a.azVx, &a.axVz, &a.azVz, &a.axSxx, &a.azSxz, &a.axSxz, &a.azSzz})
        v->assign(static_cast<size_t>(s.n_sponge_cell) + 1, 0.0f);
    a.order_x.allocate(s.ibeg, s.iend, -1);
    a.order_x(1) = a.order_x(2) = 1;
    a.order_x(s.nx - 1) = a.order_x(s.nx) = 1;
}

inline int pml_sign(const State& s, const Absorber& a, int k, int i) {
    return std::max(order_sign(s, k, i), a.order_x(i));
}

void pml_update_vel_core(State& s, Absorber& a, const Box& bb) {
    const FdCoef c = make_fd_coef(s);
    const float dt = s.dt;
#pragma omp parallel for schedule(dynamic)
    for (int i = bb.ib; i <= bb.ie; ++i) {
        for (int k = bb.kb; k <= bb.ke; ++k) {
            const int p = bb.offset + (i - bb.ib) * bb.nz + (k - bb.kb + 1);
            const int isign = pml_sign(s, a, k, i);
            const real_mp re40x = c.rc40x + isign * c.rd40x, re41x = c.rc41x + isign * c.rd41x;
            const real_mp re40z = c.rc40z + isign * c.rd40z, re41z = c.rc41z + isign * c.rd41z;
            const real_mp dxSxx = (s.Sxx(k, i + 1) - s.Sxx(k, i)) * re40x - (s.Sxx(k, i + 2) - s.Sxx(k, i - 1)) * re41x;
            const real_mp dzSzz = (s.Szz(k + 1, i) - s.Szz(k, i)) * re40z - (s.Szz(k + 2, i) - s.Szz(k - 1, i)) * re41z;
            const real_mp dxSxz = (s.Sxz(k, i) - s.Sxz(k, i - 1)) * re40x - (s.Sxz(k, i + 1) - s.Sxz(k, i - 2)) * re41x;
            const real_mp dzSxz = (s.Sxz(k, i) - s.Sxz(k - 1, i)) * re40z - (s.Sxz(k + 1, i) - s.Sxz(k - 2, i)) * re41z;
            const float bx = 2.0f / (s.rho(k, i) + s.rho(k, i + 1));
            const float bz = 2.0f / (s.rho(k, i) + s.rho(k + 1, i));
            s.Vx(k, i) = static_cast<float>(
                s.Vx(k, i) + bx * (a.gxe(1, i) * dxSxx + a.gzc(1, k) * dzSxz + a.gxe(2, i) * a.axSxx[p] +
                                   a.gzc(2, k) * a.azSxz[p]) * dt);
            s.Vz(k, i) = static_cast<float>(
                s.Vz(k, i) + bz * (a.gxc(1, i) * dxSxz + a.gze(1, k) * dzSzz + a.gxc(2, i) * a.axSxz[p] +
                                   a.gze(2, k) * a.azSzz[p]) * dt);
            a.axSxx[p] = static_cast<float>(a.gxe(3, i) * a.axSxx[p] + a.gxe(4, i) * dxSxx * dt);
            a.azSxz[p] = static_cast<float>(a.gzc(3, k) * a.azSxz[p] + a.gzc(4, k) * dzSxz * dt);
            a.axSxz[p] = static_cast<float>(a.gxc(3, i) * a.axSxz[p] + a.gxc(4, i) * dxSxz * dt);
            a.azSzz[p] = static_cast<float>(a.gze(3, k) * a.azSzz[p] + a.gze(4, k) * dzSzz * dt);
        }
    }
}

void pml_update_stress_core(State& s, Absorber& a, const Box& bb) {
    const FdCoef c = make_fd_coef(s);
    const float dt = s.dt;
    const float d2 = s.d2;
#pragma omp parallel for schedule(dynamic)
    for (int i = bb.ib; i <= bb.ie; ++i) {
        for (int k = bb.kb; k <= bb.ke; ++k) {
            const int p = bb.offset + (i - bb.ib) * bb.nz + (k - bb.kb + 1);
            const int isign = pml_sign(s, a, k, i);
            const real_mp re40x = c.rc40x + isign * c.rd40x, re41x = c.rc41x + isign * c.rd41x;
            const real_mp re40z = c.rc40z + isign * c.rd40z, re41z = c.rc41z + isign * c.rd41z;
            const real_mp dxVx = (s.Vx(k, i) - s.Vx(k, i - 1)) * re40x - (s.Vx(k, i + 1) - s.Vx(k, i - 2)) * re41x;
            const real_mp dzVz = (s.Vz(k, i) - s.Vz(k - 1, i)) * re40z - (s.Vz(k + 1, i) - s.Vz(k - 2, i)) * re41z;
            const float mu2 = 2 * s.mu(k, i);
            const float lam2mu = s.lam(k, i) + mu2;
            const float taup1 = s.taup(k, i);
            const float taus1 = s.taus(k, i);
            const float dxVx_ade = static_cast<float>(a.gxc(1, i) * dxVx + a.gxc(2, i) * a.axVx[p]);
            const float dzVz_ade = static_cast<float>(a.gzc(1, k) * dzVz + a.gzc(2, k) * a.azVz[p]);
            const float d2v2 = dxVx_ade + dzVz_ade;
            const float f_Rxx = lam2mu * taup1 * d2v2 - mu2 * taus1 * dzVz_ade;
            const float f_Rzz = lam2mu * taup1 * d2v2 - mu2 * taus1 * dxVx_ade;
            float Rxx_n = 0.0f, Rzz_n = 0.0f;
            for (int m = 1; m <= NM; ++m) {
                s.Rxx(m, k, i) = s.c1[m - 1] * s.Rxx(m, k, i) - s.c2[m - 1] * f_Rxx * dt;
                s.Rzz(m, k, i) = s.c1[m - 1] * s.Rzz(m, k, i) - s.c2[m - 1] * f_Rzz * dt;
                Rxx_n += s.d1[m - 1] * s.Rxx(m, k, i);
                Rzz_n += s.d1[m - 1] * s.Rzz(m, k, i);
            }
            const float taup_plus1 = 1 + taup1 * (1 + d2);
            const float taus_plus1 = 1 + taus1 * (1 + d2);
            s.Sxx(k, i) += (lam2mu * taup_plus1 * d2v2 - mu2 * taus_plus1 * dzVz_ade + Rxx_n) * dt;
            s.Szz(k, i) += (lam2mu * taup_plus1 * d2v2 - mu2 * taus_plus1 * dxVx_ade + Rzz_n) * dt;
            a.axVx[p] = static_cast<float>(a.gxc(3, i) * a.axVx[p] + a.gxc(4, i) * dxVx * dt);
            a.azVz[p] = static_cast<float>(a.gzc(3, k) * a.azVz[p] + a.gzc(4, k) * dzVz * dt);
        }
    }
#pragma omp parallel for schedule(dynamic)
    for (int i = bb.ib; i <= bb.ie; ++i) {
        for (int k = bb.kb; k <= bb.ke; ++k) {
            const int p = bb.offset + (i - bb.ib) * bb.nz + (k - bb.kb + 1);
            const int isign = pml_sign(s, a, k, i);
            const real_mp re40x = c.rc40x + isign * c.rd40x, re41x = c.rc41x + isign * c.rd41x;
            const real_mp re40z = c.rc40z + isign * c.rd40z, re41z = c.rc41z + isign * c.rd41z;
            const real_mp dxVz = (s.Vz(k, i + 1) - s.Vz(k, i)) * re40x - (s.Vz(k, i + 2) - s.Vz(k, i - 1)) * re41x;
            const real_mp dzVx = (s.Vx(k + 1, i) - s.Vx(k, i)) * re40z - (s.Vx(k + 2, i) - s.Vx(k - 1, i)) * re41z;
            const float taus1 = s.taus(k, i);
            const float mu_xz = harmonic_mu(s.mu(k, i), s.mu(k + 1, i), s.mu(k, i + 1), s.mu(k + 1, i + 1));
            const float dxVz_ade = static_cast<float>(a.gxe(1, i) * dxVz + a.gxe(2, i) * a.axVz[p]);
            const float dzVx_ade = static_cast<float>(a.gze(1, k) * dzVx + a.gze(2, k) * a.azVx[p]);
            const float dxVz_dzVx = dxVz_ade + dzVx_ade;
            const float f_Rxz = mu_xz * taus1 * dxVz_dzVx;
            float Rxz_n = 0.0f;
            for (int m = 1; m <= NM; ++m) {
                s.Rxz(m, k, i) = s.c1[m - 1] * s.Rxz(m, k, i) - s.c2[m - 1] * f_Rxz * dt;
                Rxz_n += s.d1[m - 1] * s.Rxz(m, k, i);
            }
            const float taus_plus1 = 1 + taus1 * (1 + d2);
            s.Sxz(k, i) += (mu_xz * taus_plus1 * dxVz_dzVx + Rxz_n) * dt;
            a.azVx[p] = static_cast<float>(a.gze(3, k) * a.azVx[p] + a.gze(4, k) * dzVx * dt);
            a.axVz[p] = static_cast<float>(a.gxe(3, i) * a.axVz[p] + a.gxe(4, i) * dxVz * dt);
        }
    }
}

void cerjan_setup(State& s, Absorber& a) {
    const float alpha = 0.09f;
    a.gx_c.allocate(s.ibeg_m, s.iend_m, 1.0f);
    a.gx_b.allocate(s.ibeg_m, s.iend_m, 1.0f);
    a.gz_c.allocate(s.kbeg_m, s.kend_m, 1.0f);
    a.gz_b.allocate(s.kbeg_m, s.kend_m, 1.0f);
    const float Lx = static_cast<float>(s.na * s.dx);
    const float Lz = static_cast<float>(s.na * s.dz);
    const float dx = static_cast<float>(s.dx), dz = static_cast<float>(s.dz);
    auto prof = [&](float v, float L) {
        const float r = 1.0f - v / L;
        return std::exp(-alpha * (r * r));  // Fortran: -alpha * (...)**2
    };
    for (int i = s.ibeg; i <= s.iend; ++i) {
        if (i <= s.na) {
            a.gx_c(i) = prof(i2x(i, 0.0f, dx), Lx);
            a.gx_b(i) = prof(i2x(i, 0.0f, dx) + dx / 2, Lx);
        } else if (i >= s.nx - s.na + 1) {
            a.gx_c(i) = prof(i2x(i, s.nx * dx, -dx) + dx / 2, Lx);
            a.gx_b(i) = prof(i2x(i, s.nx * dx, -dx), Lx);
        }
    }
    for (int k = s.kbeg; k <= s.kend; ++k) {
        if (k <= s.na) {
            if (s.fullspace_mode) {
                a.gz_c(k) = prof(k2z(k, 0.0f, dz), Lz);
                a.gz_b(k) = prof(k2z(k, 0.0f, dz) + dz / 2, Lz);
            }
        } else if (k >= s.nz - s.na + 1) {
            a.gz_c(k) = prof(k2z(k, s.nz * dz, -dz) + dz / 2, Lz);
            a.gz_b(k) = prof(k2z(k, s.nz * dz, -dz), Lz);
        }
    }
}

}  // namespace

void absorb_setup(State& s, Absorber& a) {
    if (s.abc_type == "pml") pml_setup(s, a);
    else if (s.abc_type == "cerjan") cerjan_setup(s, a);
    else throw std::runtime_error("abc_type must be 'pml' or 'cerjan': " + s.abc_type);
}

void absorb_update_stress(State& s, Absorber& a) {
    if (s.abc_type == "pml") {
        for (const Box& b : s.box)
            if (b.ncell > 0) pml_update_stress_core(s, a, b);
    } else {
#pragma omp parallel for schedule(dynamic)
        for (int i = s.ibeg; i <= s.iend; ++i)
            for (int k = s.kbeg; k <= s.kend; ++k) {
                const float gc = a.gx_c(i) * a.gz_c(k);
                const float gb = a.gx_b(i) * a.gz_b(k);
                s.Sxx(k, i) *= gc;
                s.Szz(k, i) *= gc;
                s.Sxz(k, i) *= gb;
            }
    }
}

void absorb_update_vel(State& s, Absorber& a) {
    if (s.abc_type == "pml") {
        for (const Box& b : s.box)
            if (b.ncell > 0) pml_update_vel_core(s, a, b);
    } else {
#pragma omp parallel for schedule(dynamic)
        for (int i = s.ibeg; i <= s.iend; ++i)
            for (int k = s.kbeg; k <= s.kend; ++k) {
                s.Vx(k, i) = s.Vx(k, i) * a.gx_b(i) * a.gz_c(k);
                s.Vz(k, i) = s.Vz(k, i) * a.gx_c(i) * a.gz_b(k);
            }
    }
}

}  // namespace swpc
