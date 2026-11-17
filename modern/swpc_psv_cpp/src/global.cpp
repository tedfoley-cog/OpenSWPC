#include <algorithm>
#include <cmath>
#include <ctime>
#include <filesystem>

#include "fdtool.hpp"
#include "global.hpp"
#include "readini.hpp"
#include "state.hpp"

namespace swpc {

namespace {
void set_box(Box& b, int ib, int ie, int kb, int ke, int& offset) {
    b.ib = ib; b.ie = ie; b.kb = kb; b.ke = ke;
    b.nx = std::max(ie - ib + 1, 0);
    b.nz = std::max(ke - kb + 1, 0);
    b.ncell = b.nx * b.nz;
    b.offset = offset;
    offset += b.ncell;
}
}  // namespace

void global_setup(State& s, const IniFile& prm) {
    s.benchmark_mode = prm.get_bool("benchmark_mode", false);
    s.title = prm.get_string("title", "swpc_psv");
    s.nproc_x = prm.get_int("nproc_x", 1);
    s.nx = prm.get_int("nx", 256);
    s.nz = prm.get_int("nz", 256);
    s.nt = prm.get_int("nt", 1000);
    s.ipad = prm.get_int("ipad", 0);
    s.kpad = prm.get_int("kpad", 0);
    s.odir = prm.get_string("odir", "./out");
    if (s.benchmark_mode) {
        s.dx = 0.5; s.dz = 0.5; s.dt = 0.04f; s.na = 20;
        s.xbeg = -s.nx / 2.0f * static_cast<float>(s.dx);
        s.zbeg = -30 * static_cast<float>(s.dz);
        s.tbeg = 0.0f; s.clon = 139.7604f; s.clat = 35.7182f; s.phi = 0.0f;
        s.abc_type = "pml"; s.fullspace_mode = false;
    } else {
        s.dx = prm.get_double("dx", 0.5);
        s.dz = prm.get_double("dz", 0.5);
        s.dt = prm.get_float("dt", 0.01f);
        s.na = prm.get_int("na", 20);
        s.xbeg = prm.get_float("xbeg", -(s.nx / 2) * static_cast<float>(s.dx));
        s.zbeg = prm.get_float("zbeg", -30 * static_cast<float>(s.dz));
        s.tbeg = prm.get_float("tbeg", 0.0f);
        s.clon = prm.get_float("clon", 139.7604f);
        s.clat = prm.get_float("clat", 35.7182f);
        s.phi = prm.get_float("phi", 0.0f);
        s.abc_type = prm.get_string("abc_type", "pml");
        s.fullspace_mode = prm.get_bool("fullspace_mode", false);
    }
    s.exedate = static_cast<int>(std::time(nullptr));
    s.xend = s.xbeg + s.nx * static_cast<float>(s.dx);
    s.zend = s.zbeg + s.nz * static_cast<float>(s.dz);
}

void global_setup2(State& s) {
    // single rank: nproc_x is accepted for input compatibility but the whole
    // domain i = 1..nx is owned by this process.
    s.myid = 0;
    s.idx = 0;
    s.nproc_x = 1;
    std::error_code ec;
    std::filesystem::create_directories(s.odir, ec);

    s.ibeg = 1;
    s.iend = s.nx;
    s.kbeg = 1;
    s.kend = s.nz;
    s.ibeg_m = s.ibeg - 3;
    s.iend_m = s.iend + 3 + s.ipad;
    s.kbeg_m = s.kbeg - 3;
    s.kend_m = s.kend + 3 + s.kpad;

    s.xc.allocate(s.ibeg_m, s.iend_m);
    s.zc.allocate(s.kbeg_m, s.kend_m);
    for (int i = s.ibeg_m; i <= s.iend_m; ++i) s.xc(i) = i2x(i, s.xbeg, static_cast<float>(s.dx));
    for (int k = s.kbeg_m; k <= s.kend_m; ++k) s.zc(k) = k2z(k, s.zbeg, static_cast<float>(s.dz));

    const int na = s.na, nx = s.nx;
    s.ibeg_k = s.ibeg; s.iend_k = s.iend; s.kbeg_k = s.kbeg; s.kend_k = s.kend;
    if (s.abc_type == "pml") {
        if (s.fullspace_mode) s.kbeg_k = na + 1;
        if (s.iend <= na) s.ibeg_k = s.iend + 1;
        else if (s.ibeg <= na) s.ibeg_k = na + 1;
        if (s.ibeg >= nx - na + 1) s.iend_k = s.ibeg - 1;
        else if (s.iend >= nx - na + 1) s.iend_k = nx - na;
        s.kend_k = s.nz - na;
    }
    int offset = 1;
    const int ib = s.ibeg, ie = s.iend, kb = s.kbeg, ke = s.kend;
    set_box(s.box[0], ib, na, kb, ke, offset);
    set_box(s.box[1], nx - na + 1, ie, kb, ke, offset);
    if (s.fullspace_mode)
        set_box(s.box[2], std::max(na + 1, ib), std::min(nx - na, ie), kb, na, offset);
    else
        set_box(s.box[2], std::max(na + 1, ib), std::min(nx - na, ie), kb, kb - 1, offset);
    set_box(s.box[3], std::max(na + 1, ib), std::min(nx - na, ie), ke - na + 1, ke, offset);
    s.n_sponge_cell = offset;
}

}  // namespace swpc
