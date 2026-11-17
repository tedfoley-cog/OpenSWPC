// Port of src/shared/m_fdtool.f90 (routines used by swpc_psv).
// Fortran real(SP) -> float; PI is double precision (m_std), so mixed
// expressions are evaluated in double and rounded to float exactly as gfortran does.
#pragma once
#include <string>
#include <vector>

namespace swpc {

constexpr double PI = 3.141592653589793238462643383279502884;

float moment_magnitude(float m0);
float seismic_moment(float mw);
void sdr2moment(float strike, float dip, float rake, float& mxx, float& myy, float& mzz,
                float& myz, float& mxz, float& mxy);

float kupper(float t, float ts, float tr);
float texp(float t, float ts, float tr);
float cosine(float t, float ts, float tr);
float boxcar(float t, float ts, float tr);
float triangle(float t, float ts, float tr);
float herrmann(float t, float ts, float tr);
float asymcos(float t, float ts, float tr1, float tr2);
float momentrate(float t, const std::string& stftype, const float* srcprm);

int x2i(float x, float xbeg, float dx);
int z2k(float z, float zbeg, float dz);
float i2x(int i, float xbeg, float dx);
float k2z(int k, float zbeg, float dz);

void visco_set_relaxtime(int nm, std::vector<float>& ts, float fmin, float fmax);
float visco_chi(int nm, const std::vector<float>& ts, float tau, float f);
float visco_constq_zeta(int nm, float fmin, float fmax, const std::vector<float>& ts);

void fdm_stable_dt(float dx, float dy, float dz, float vmax, float& dt);
void fdm_cond_stability(float dx, float dy, float dz, float vmax, float dt, float& c);
void fdm_cond_wavelength(float dx, float dy, float dz, float vmin, float fmax, float& r);
float deg2rad(float deg);
double deg2rad(double deg);

}  // namespace swpc
