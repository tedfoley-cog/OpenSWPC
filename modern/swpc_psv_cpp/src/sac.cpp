#include "sac.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace swpc {

SacHeader::SacHeader() {
    for (float& v : f) v = -12345.0f;
    for (int& v : i) v = -12345;
    for (int& v : l) v = 0;
    for (auto& s : k) s = "-12345";
    i[NVHDR - 71] = 6;
    i[IFTYPE - 71] = 1;
    l[LEVEN - 106] = 1;
    l[LOVROK - 106] = 1;
}

namespace {
std::string pad(const std::string& s, size_t n) {
    std::string r = s.substr(0, std::min(s.size(), n));
    r.resize(n, ' ');
    return r;
}
}  // namespace

void sac_write(const std::string& fn, const SacHeader& sh, const std::vector<float>& dat) {
    float fh[70];
    std::copy(sh.f, sh.f + 70, fh);
    fh[0] = static_cast<float>(static_cast<long long>(sh.f[0] * 1e7)) / 1e7f;

    std::ofstream out(fn, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("cannot write SAC file: " + fn);
    out.write(reinterpret_cast<const char*>(fh), sizeof(fh));
    out.write(reinterpret_cast<const char*>(sh.i), sizeof(sh.i));
    out.write(reinterpret_cast<const char*>(sh.l), sizeof(sh.l));
    for (int j = 0; j < 23; ++j) out << pad(sh.k[j], j == KEVNM ? 16 : 8);
    out.write(reinterpret_cast<const char*>(dat.data()), static_cast<std::streamsize>(dat.size() * sizeof(float)));
    out.close();
    if (!out) throw std::runtime_error("failed writing SAC file: " + fn);
}

}  // namespace swpc
