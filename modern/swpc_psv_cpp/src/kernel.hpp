// Port of m_kernel.f90 (interior 4th-order staggered-grid FDM) and
// m_absorb.f90 / m_absorb_p.f90 / m_absorb_c.f90 (PML and Cerjan absorbers).
#pragma once
#include "readini.hpp"
#include "state.hpp"

namespace swpc {

struct FdCoef {
    real_mp rc40x, rc40z, rc41x, rc41z, rd40x, rd40z, rd41x, rd41z;
};
FdCoef make_fd_coef(const State& s);

void kernel_setup(State& s);
void kernel_update_vel(State& s);
void kernel_update_stress(State& s);

struct Absorber {
    // PML
    Array2D<float> gxc, gxe, gzc, gze;  // (1:4, i) / (1:4, k)
    std::vector<float> axVx, azVx, axVz, azVz, axSxx, azSxz, axSxz, azSzz;
    Array1D<int> order_x;
    // Cerjan
    Array1D<float> gx_c, gx_b, gz_c, gz_b;
};

void absorb_setup(State& s, Absorber& a);
void absorb_update_stress(State& s, Absorber& a);
void absorb_update_vel(State& s, Absorber& a);
void damping_profile(const State& s, float x, float H, float xb, float xe, float g[4]);

}  // namespace swpc
