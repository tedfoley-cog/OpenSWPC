// Solver-wide state: port of src/swpc_psv/m_global.f90 (single-node; the MPI
// decomposition collapses to one rank owning i = 1..nx).
#pragma once
#include <algorithm>
#include <array>
#include <string>
#include <vector>

namespace swpc {

// Fortran-ordered 2D array a(kb:ke, ib:ie); k is the contiguous index.
template <class T>
class Array2D {
public:
    Array2D() = default;
    void allocate(int kb, int ke, int ib, int ie, T v = T()) {
        kb_ = kb; ke_ = ke; ib_ = ib; ie_ = ie;
        nk_ = ke - kb + 1;
        data_.assign(static_cast<size_t>(nk_) * static_cast<size_t>(ie - ib + 1), v);
    }
    T& operator()(int k, int i) { return data_[static_cast<size_t>(k - kb_) + static_cast<size_t>(i - ib_) * nk_]; }
    const T& operator()(int k, int i) const {
        return data_[static_cast<size_t>(k - kb_) + static_cast<size_t>(i - ib_) * nk_];
    }
    bool empty() const { return data_.empty(); }
    int kb() const { return kb_; }
    int ke() const { return ke_; }
    int ib() const { return ib_; }
    int ie() const { return ie_; }
    std::vector<T>& raw() { return data_; }
    const std::vector<T>& raw() const { return data_; }

private:
    int kb_ = 0, ke_ = -1, ib_ = 0, ie_ = -1, nk_ = 0;
    std::vector<T> data_;
};

// Fortran-ordered a(1:nm, kb:ke, ib:ie) for viscoelastic memory variables.
template <class T>
class Array3D {
public:
    void allocate(int nm, int kb, int ke, int ib, int ie, T v = T()) {
        nm_ = nm; kb_ = kb; ib_ = ib; nk_ = ke - kb + 1;
        data_.assign(static_cast<size_t>(nm) * nk_ * static_cast<size_t>(ie - ib + 1), v);
    }
    T& operator()(int m, int k, int i) {
        return data_[static_cast<size_t>(m - 1) + nm_ * (static_cast<size_t>(k - kb_) + static_cast<size_t>(i - ib_) * nk_)];
    }

private:
    int nm_ = 0, kb_ = 0, ib_ = 0, nk_ = 0;
    std::vector<T> data_;
};

// 1D array with arbitrary lower bound, a(lb:ub).
template <class T>
class Array1D {
public:
    void allocate(int lb, int ub, T v = T()) { lb_ = lb; data_.assign(static_cast<size_t>(ub - lb + 1), v); }
    T& operator()(int i) { return data_[static_cast<size_t>(i - lb_)]; }
    const T& operator()(int i) const { return data_[static_cast<size_t>(i - lb_)]; }
    void fill(T v) { std::fill(data_.begin(), data_.end(), v); }
    int lb() const { return lb_; }
    int ub() const { return lb_ + static_cast<int>(data_.size()) - 1; }

private:
    int lb_ = 0;
    std::vector<T> data_;
};

struct Box {
    int ib = 0, ie = -1, nx = 0;
    int kb = 0, ke = -1, nz = 0;
    int ncell = 0;
    int offset = 0;
};

using real_mp = double;  // Fortran MP (mixed precision) == DP

constexpr int NM = 3;    // number of viscoelastic memory variables
constexpr int NBD = 9;   // number of memorised boundary depths

struct State {
    float UC = 1.0e-12f;  // conventional -> SI unit for moment tensor

    Array2D<float> Vx, Vz;
    Array2D<real_mp> Sxx, Szz, Sxz;
    Array3D<float> Rxx, Rzz, Rxz;
    Array2D<float> rho, lam, mu, taup, taus;
    std::vector<float> ts;

    std::string title = "swpc_psv";
    int exedate = 0;
    int nx = 256, nz = 256, nt = 1000;
    real_mp dx = 0.5, dz = 0.5;
    float dt = 0.01f;
    float xbeg = 0, xend = 0, zbeg = 0, zend = 0, tbeg = 0, tend = 0;
    float vmin = 0, vmax = 0, fmax = 0, fcut = 0;
    int nproc_x = 1, myid = 0, idx = 0;
    int ibeg = 1, iend = 1, kbeg = 1, kend = 1;
    int ibeg_m = 0, iend_m = 0, kbeg_m = 0, kend_m = 0;
    int ipad = 0, kpad = 0;
    float M0 = 0;
    int na = 20;
    int ibeg_k = 0, iend_k = 0, kbeg_k = 0, kend_k = 0;
    std::string abc_type = "pml";
    std::string odir = "./out";
    Array1D<int> kfs, kob, kfs_top, kfs_bot, kob_top, kob_bot;
    Array2D<float> bddep;  // bddep(i, 0:NBD) stored as (id, i)
    float clon = 139.7604f, clat = 35.7182f, phi = 0.0f;
    Array1D<float> xc, zc;
    float evlo = 0, evla = 0, evdp = 0;
    float mxx0 = 0, myy0 = 0, mzz0 = 0, myz0 = 0, mxz0 = 0, mxy0 = 0;
    float fx0 = 0, fy0 = 0, fz0 = 0;
    float otim = 0;
    float sx0 = 0, sy0 = 0;
    bool fullspace_mode = false, benchmark_mode = false, pw_mode = false, bf_mode = false;
    std::array<Box, 4> box{};
    int n_sponge_cell = 0;
    std::vector<float> c1, c2, d1;
    float d2 = 0;
};

class IniFile;

void global_setup(State& s, const IniFile& prm);
void global_setup2(State& s);

}  // namespace swpc
