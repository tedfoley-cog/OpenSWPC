#include "fdtool.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <complex>

namespace swpc {

float deg2rad(float deg) { return static_cast<float>(PI / 180.0 * static_cast<double>(deg)); }
double deg2rad(double deg) { return PI / 180.0 * deg; }

float moment_magnitude(float m0) {
    if (m0 < FLT_EPSILON) return -12345.0f;
    return (std::log10(m0) - 9.1f) * 2.0f / 3.0f;
}

float seismic_moment(float mw) { return std::pow(10.0f, 1.5f * mw + 9.05f); }

void sdr2moment(float strike, float dip, float rake, float& mxx, float& myy, float& mzz,
                float& myz, float& mxz, float& mxy) {
    const float sind = std::sin(deg2rad(dip));
    const float cosd = std::cos(deg2rad(dip));
    const float sin2d = std::sin(deg2rad(2 * dip));
    const float cos2d = std::cos(deg2rad(2 * dip));
    const float sinl = std::sin(deg2rad(rake));
    const float cosl = std::cos(deg2rad(rake));
    const float sinf = std::sin(deg2rad(strike));
    const float cosf = std::cos(deg2rad(strike));
    const float sin2f = std::sin(deg2rad(2 * strike));
    const float cos2f = std::cos(deg2rad(2 * strike));
    mxx = -(sind * cosl * sin2f + sin2d * sinl * sinf * sinf);
    mxy = (sind * cosl * cos2f + sin2d * sinl * sin2f / 2);
    mxz = -(cosd * cosl * cosf + cos2d * sinl * sinf);
    myy = (sind * cosl * sin2f - sin2d * sinl * cosf * cosf);
    myz = -(cosd * cosl * sinf - cos2d * sinl * cosf);
    mzz = (sin2d * sinl);
}

float kupper(float t, float ts, float tr) {
    if (ts <= t && t <= ts + tr) {
        const double s = std::sin(PI * (t - ts) / tr);
        return static_cast<float>(3 * PI * (s * s * s) / static_cast<double>(4 * tr));
    }
    return 0.0f;
}

float texp(float t, float ts, float tr) {
    if (ts <= t) {
        const float tt = t - ts;
        return static_cast<float>((2 * PI) * (2 * PI) * tt / static_cast<double>(tr * tr) *
                                  std::exp(-2 * PI * tt / tr));
    }
    return 0.0f;
}

float cosine(float t, float ts, float tr) {
    if (ts <= t && t <= ts + tr) return static_cast<float>((1 - std::cos(2 * PI * (t - ts) / tr)) / tr);
    return 0.0f;
}

float boxcar(float t, float ts, float tr) {
    if (ts <= t && t <= ts + tr) return 1.0f / tr;
    return 0.0f;
}

float triangle(float t, float ts, float tr) {
    if (ts <= t && t <= ts + tr / 2) return 4 * (t - ts) / (tr * tr);
    if (ts + tr / 2 < t && t <= ts + tr) return -4 * (t - ts - tr) / (tr * tr);
    return 0.0f;
}

float herrmann(float t, float ts, float tr) {
    const float t1 = ts + tr / 4;
    const float t2 = ts + 3 * tr / 4;
    const float tr3 = tr * tr * tr;
    if (ts <= t && t < t1) return 16 * (t - ts) * (t - ts) / tr3;
    if (t1 <= t && t < t2)
        return -2 * (8 * (t * t + tr * ts + ts * ts - t * tr - 2 * t * ts) + tr * tr) / tr3;
    if (t2 <= t && t <= ts + tr) return 16 * (ts + tr - t) * (ts + tr - t) / tr3;
    return 0.0f;
}

float asymcos(float t, float ts, float tr1, float tr2) {
    const float tr = tr1 + tr2;
    if (ts <= t && t <= ts + tr1) return static_cast<float>((1.0 - std::cos(PI * (t - ts) / tr1)) / tr);
    if (ts + tr1 < t && t <= ts + tr1 + tr2)
        return static_cast<float>((1.0 + std::cos(PI * (t - ts - tr1) / tr2)) / tr);
    return 0.0f;
}

float momentrate(float t, const std::string& stftype, const float* p) {
    if (stftype == "boxcar") return boxcar(t, p[0], p[1]);
    if (stftype == "triangle") return triangle(t, p[0], p[1]);
    if (stftype == "herrmann") return herrmann(t, p[0], p[1]);
    if (stftype == "kupper") return kupper(t, p[0], p[1]);
    if (stftype == "cosine") return cosine(t, p[0], p[1]);
    if (stftype == "texp") return texp(t, p[0], p[1]);
    if (stftype == "asymcos") return asymcos(t, p[0], p[1], p[2]);
    return kupper(t, p[0], p[1]);
}

int x2i(float x, float xbeg, float dx) { return static_cast<int>(std::ceil((x - xbeg) / dx)); }
int z2k(float z, float zbeg, float dz) { return static_cast<int>(std::ceil((z - zbeg) / dz)); }
float i2x(int i, float xbeg, float dx) { return xbeg + (static_cast<float>(i) - 0.5f) * dx; }
float k2z(int k, float zbeg, float dz) { return zbeg + (static_cast<float>(k) - 0.5f) * dz; }

void visco_set_relaxtime(int nm, std::vector<float>& ts, float fmin, float fmax) {
    const float omega_a = static_cast<float>(2 * PI * fmin);
    const float omega_b = static_cast<float>(2 * PI * fmax);
    if (nm == 0) return;
    if (nm == 1) {
        const float omega = std::sqrt(omega_a * omega_b);
        ts[0] = 1.0f / omega;
        return;
    }
    for (int im = 1; im <= nm; ++im) {
        const float omega = static_cast<float>(
            omega_a * std::pow(static_cast<double>(omega_b / omega_a), double(im - 1) / double(nm - 1)));
        ts[im - 1] = 1.0f / omega;
    }
}

float visco_chi(int nm, const std::vector<float>& ts, float tau, float f) {
    const std::complex<double> EI(0.0, 1.0);
    const float omega = static_cast<float>(2 * PI * f);
    std::complex<float> cc(0.0f, 0.0f);
    for (int im = 0; im < nm; ++im) {
        const std::complex<double> num = 1.0 - EI * static_cast<double>(omega) * static_cast<double>(ts[im]) *
                                                   static_cast<double>(1 + tau);
        const std::complex<double> den = 1.0 - EI * static_cast<double>(omega) * static_cast<double>(ts[im]);
        cc = std::complex<float>(std::complex<double>(cc) + num / den);
    }
    cc = cc / static_cast<float>(nm);
    return 1.0f / std::real(std::pow(cc, -0.5f));
}

float visco_constq_zeta(int nm, float fmin, float fmax, const std::vector<float>& ts) {
    if (nm == 0) return 0.0f;
    const float om_a = static_cast<float>(2 * PI * fmin);
    const float om_b = static_cast<float>(2 * PI * fmax);
    float i0sum = 0.0f, i1sum = 0.0f, i2sum = 0.0f;
    for (int im = 0; im < nm; ++im) {
        const float t = ts[im];
        i0sum += (std::log(1.0f + om_b * om_b * t * t) - std::log(1.0f + om_a * om_a * t * t)) / (2 * t);
        i1sum += ((std::atan(om_b * t) - om_b * t / (1 + om_b * om_b * t * t)) -
                  (std::atan(om_a * t) - om_a * t / (1 + om_a * om_a * t * t))) /
                 (2 * t);
    }
    for (int im = 0; im < nm - 1; ++im) {
        for (int km = im + 1; km < nm; ++km) {
            const float wk1 = std::atan(om_b * ts[im]) / ts[im] - std::atan(om_b * ts[km]) / ts[km];
            const float wk2 = std::atan(om_a * ts[im]) / ts[im] - std::atan(om_a * ts[km]) / ts[km];
            i2sum += ts[im] * ts[km] / (ts[km] * ts[km] - ts[im] * ts[im]) * (wk1 - wk2);
        }
    }
    return i0sum / (i1sum + 2 * i2sum);
}

void fdm_stable_dt(float dx, float dy, float dz, float vmax, float& dt) {
    const float hh = static_cast<float>(1. / std::sqrt(1 / (dx * dx) + 1 / (dy * dy) + 1 / (dz * dz)));
    const float cc = static_cast<float>(6. / 7.);
    dt = cc * hh / vmax;
}

void fdm_cond_stability(float dx, float dy, float dz, float vmax, float dt, float& c) {
    float dt2;
    fdm_stable_dt(dx, dy, dz, vmax, dt2);
    c = dt / dt2;
}

void fdm_cond_wavelength(float dx, float dy, float dz, float vmin, float fmax, float& r) {
    const float dh = std::max({dx, dy, dz});
    const float lambda_min = vmin / fmax;
    r = lambda_min / dh;
}

}  // namespace swpc
