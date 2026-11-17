#include "geomap.hpp"

#include <cmath>

#include "fdtool.hpp"

namespace swpc {

namespace {
constexpr double a = 6378137.0;
constexpr double F = 298.257222101;
constexpr double m0 = 0.9999;
constexpr double n = 1.0 / (2.0 * F - 1.0);

struct Coef {
    double alpha[6], beta[6], AA[6], delta[7];
    Coef() {
        alpha[1] = (1 / 2. + (-2 / 3. + (5 / 16. + (41 / 180. - 127 / 288. * n) * n) * n) * n) * n;
        alpha[2] = (13 / 48. + (-3 / 5. + (557 / 1440. + 281 / 630. * n) * n) * n) * std::pow(n, 2);
        alpha[3] = (61 / 240. + (-103 / 140. + 15061 / 26880. * n) * n) * std::pow(n, 3);
        alpha[4] = (49561 / 161280. - 179 / 168. * n) * std::pow(n, 4);
        alpha[5] = 34729 / 80640. * std::pow(n, 5);
        beta[1] = (1 / 2. + (-2 / 3. + (37 / 96. + (-1 / 360. - 81 / 512. * n) * n) * n) * n) * n;
        beta[2] = ((1 / 48. + (1 / 15. + (-437 / 1440. + 46 / 105. * n) * n) * n) * n) * n;
        beta[3] = (((17 / 480. + (-37 / 840. - 209 / 4480. * n) * n) * n) * n) * n;
        beta[4] = ((((4397 / 161280. - 11 / 504. * n) * n) * n) * n) * n;
        beta[5] = ((((4583 / 161280. * n) * n) * n) * n) * n;
        delta[1] = (2 / 1. + (-2 / 3. + (-2 / 1. + (116 / 45. + (26 / 45. + (-2854 / 675.) * n) * n) * n) * n) * n) * n;
        delta[2] = ((7 / 3. + (-8 / 5. + (-227 / 45. + (2704 / 315. + (2323 / 945.) * n) * n) * n) * n) * n) * n;
        delta[3] = (((56 / 15. + (-136 / 35. + (-1262 / 105. + (73814 / 2835.) * n) * n) * n) * n) * n) * n;
        delta[4] = ((((4279 / 630. + (-332 / 35. + (-399572 / 14175.) * n) * n) * n) * n) * n) * n;
        delta[5] = (((((4174 / 315. + (-144838 / 6237.) * n) * n) * n) * n) * n) * n;
        delta[6] = ((((((601676 / 22275.) * n) * n) * n) * n) * n) * n;
        AA[0] = 1 + (1 / 4. + 1 / 64. * n * n) * n * n;
        AA[1] = -3 / 2. * (1. - 1 / 8. * n * n - 1 / 64. * std::pow(n, 4)) * n;
        AA[2] = 15 / 16. * (1. - 1 / 4. * n * n) * n * n;
        AA[3] = -35 / 48. * (1. - 5 / 16. * n * n) * std::pow(n, 3);
        AA[4] = 315 / 512. * std::pow(n, 4);
        AA[5] = -693 / 1280. * std::pow(n, 5);
    }
};
const Coef& coef() {
    static const Coef c;
    return c;
}

double S_phi0(double phi0) {
    const Coef& c = coef();
    double s = c.AA[0] * phi0;
    for (int j = 1; j <= 5; ++j) s += c.AA[j] * std::sin(2 * j * phi0);
    return s * (m0 * a / (1 + n));
}

void lltoxy(double lon, double lat, double lon0, double lat0, double& x, double& y) {
    const Coef& c = coef();
    const double lam = deg2rad(lon), lam0 = deg2rad(lon0), phi = deg2rad(lat), phi0 = deg2rad(lat0);
    const double e2n = 2.0 * std::sqrt(n) / (1.0 + n);
    const double lam_c = std::cos(lam - lam0), lam_s = std::sin(lam - lam0);
    const double tan_chi = std::sinh(std::atanh(std::sin(phi)) - e2n * std::atanh(e2n * std::sin(phi)));
    const double cos_chi = std::sqrt(1 + tan_chi * tan_chi);
    const double xi = std::atan(tan_chi / lam_c);
    const double eta = std::atanh(lam_s / cos_chi);
    const double Abar = m0 * a / (1 + n) * c.AA[0];
    x = xi;
    y = eta;
    for (int j = 1; j <= 5; ++j) {
        x += c.alpha[j] * std::sin(2 * j * xi) * std::cosh(2 * j * eta);
        y += c.alpha[j] * std::cos(2 * j * xi) * std::sinh(2 * j * eta);
    }
    x = (Abar * x - S_phi0(phi0)) / 1000;
    y = (Abar * y) / 1000;
}

void xytoll(double x, double y, double lon0, double lat0, double& lon, double& lat) {
    const Coef& c = coef();
    const double lam0 = deg2rad(lon0), phi0 = deg2rad(lat0);
    const double Abar = m0 * a / (1 + n) * c.AA[0];
    const double xi = (x * 1000 + S_phi0(phi0)) / Abar;
    const double eta = y * 1000 / Abar;
    double xi2 = xi, eta2 = eta;
    for (int j = 1; j <= 5; ++j) {
        xi2 -= c.beta[j] * std::sin(2 * j * xi) * std::cosh(2 * j * eta);
        eta2 -= c.beta[j] * std::cos(2 * j * xi) * std::sinh(2 * j * eta);
    }
    const double chi = std::asin(std::sin(xi2) / std::cosh(eta2));
    const double lam = lam0 + std::atan(std::sinh(eta2) / std::cos(xi2));
    double phi = chi;
    for (int j = 1; j <= 6; ++j) phi += c.delta[j] * std::sin(2 * j * chi);
    lon = 180.0 / PI * lam;
    lat = 180.0 / PI * phi;
}
}  // namespace

void geomap_g2c(float lon, float lat, float lon0, float lat0, float phi, float& x, float& y) {
    const float phi_r = deg2rad(phi);
    double xd, yd;
    lltoxy(lon, lat, lon0, lat0, xd, yd);
    const float xx = static_cast<float>(xd), yy = static_cast<float>(yd);
    x = std::cos(phi_r) * xx + std::sin(phi_r) * yy;
    y = -std::sin(phi_r) * xx + std::cos(phi_r) * yy;
}

void geomap_c2g(float x, float y, float lon0, float lat0, float phi, float& lon, float& lat) {
    const float phi_r = deg2rad(phi);
    const float xx = std::cos(phi_r) * x - std::sin(phi_r) * y;
    const float yy = std::sin(phi_r) * x + std::cos(phi_r) * y;
    double lo, la;
    xytoll(xx, yy, lon0, lat0, lo, la);
    lon = static_cast<float>(lo);
    lat = static_cast<float>(la);
}

}  // namespace swpc
