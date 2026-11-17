#include "dump.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

namespace swpc {

namespace {
const char* dump_dir() {
    const char* d = std::getenv("SWPC_DUMP_DIR");
    return (d && *d) ? d : nullptr;
}

template <class T>
void put(std::ofstream& o, const Array2D<T>& a) {
    o.write(reinterpret_cast<const char*>(a.raw().data()), static_cast<std::streamsize>(a.raw().size() * sizeof(T)));
}

void header(std::ofstream& o, const State& s) {
    const int h[4] = {s.kbeg_m, s.kend_m, s.ibeg_m, s.iend_m};
    o.write(reinterpret_cast<const char*>(h), sizeof(h));
}
}  // namespace

void dump_medium(const State& s) {
    const char* d = dump_dir();
    if (!d) return;
    std::ofstream o(std::string(d) + "/medium.bin", std::ios::binary | std::ios::trunc);
    header(o, s);
    put(o, s.rho); put(o, s.lam); put(o, s.mu); put(o, s.taup); put(o, s.taus);
}

void dump_state(const State& s, int it) {
    const char* d = dump_dir();
    if (!d) return;
    const char* e = std::getenv("SWPC_DUMP_EVERY");
    const int every = (e && *e) ? std::atoi(e) : 1;
    if (every <= 0 || it % every != 0) return;
    char name[32];
    std::snprintf(name, sizeof(name), "/state_%08d.bin", it);
    std::ofstream o(std::string(d) + name, std::ios::binary | std::ios::trunc);
    header(o, s);
    put(o, s.Vx); put(o, s.Vz); put(o, s.Sxx); put(o, s.Szz); put(o, s.Sxz);
}

}  // namespace swpc
