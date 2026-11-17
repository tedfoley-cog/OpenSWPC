#include "source.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "fdtool.hpp"
#include "geomap.hpp"

namespace swpc {

namespace {
constexpr double R_EARTH = 6371.0;

void require(bool cond, const std::string& msg) {
    if (!cond) throw std::runtime_error("source: " + msg);
}

std::vector<std::string> data_lines(const std::string& fn) {
    std::ifstream in(fn);
    require(static_cast<bool>(in), "cannot open source file '" + fn + "'");
    std::vector<std::string> out;
    std::string line;
    while (std::getline(in, line)) {
        size_t a = line.find_first_not_of(" \t\r");
        if (a == std::string::npos) continue;
        if (line[a] == '#') continue;
        out.push_back(line.substr(a));
    }
    return out;
}

// Fortran list-directed read of n reals (comma/space separated).
std::vector<double> read_reals(const std::string& line, size_t n) {
    std::string t = line;
    for (auto& c : t)
        if (c == ',') c = ' ';
    std::istringstream ss(t);
    std::vector<double> v;
    std::string tok;
    while (v.size() < n && ss >> tok) {
        for (auto& c : tok)
            if (c == 'd' || c == 'D') c = 'e';
        v.push_back(std::stod(tok));
    }
    require(v.size() == n, "too few columns in line: " + line);
    return v;
}

void check_sdr(float strike, float dip, float rake) {
    require(-360.f <= strike && strike <= 360.f, "strike out of range");
    require(-90.f <= dip && dip <= 90.f, "dip out of range");
    require(-180.f <= rake && rake <= 180.f, "rake out of range");
}

void check_ll(float lon, float lat) {
    require(-360.f <= lon && lon <= 360.f, "longitude out of range");
    require(-90.f <= lat && lat <= 90.f, "latitude out of range");
}

// Fortran real(SP)**integer: binary exponentiation in single precision (libgcc __powisf2).
float powi_sp(float x, int m) {
    unsigned n = m < 0 ? 0u - static_cast<unsigned>(m) : static_cast<unsigned>(m);
    float y = (n % 2) ? x : 1.0f;
    while (n >>= 1) {
        x = x * x;
        if (n % 2) y = y * x;
    }
    return m < 0 ? 1.0f / y : y;
}

struct RawSource {
    float sx = 0, sy = 0, sz = 0, mo = 0, mxx = 0, mzz = 0, mxz = 0, fx = 0, fz = 0;
    std::vector<float> prm;
};

void grid_moment(State& s, const std::string& fn, const std::string& fmt, int np, bool ef,
                 std::vector<RawSource>& out) {
    const auto lines = data_lines(fn);
    for (size_t n = 0; n < lines.size(); ++n) {
        const std::string& line = lines[n];
        RawSource r;
        r.prm.assign(np, 0.0f);
        auto take_prm = [&](const std::vector<double>& v, size_t off) {
            for (int j = 0; j < np; ++j) r.prm[j] = static_cast<float>(v[off + j]);
        };
        float rd1, rd2, rd3;
        const bool ll = fmt.rfind("ll", 0) == 0;
        if (fmt == "xym0ij" || fmt == "llm0ij" || fmt == "xymwij" || fmt == "llmwij") {
            auto v = read_reals(line, 3 + np + 7);
            take_prm(v, 3);
            const float m = static_cast<float>(v[3 + np]);
            r.mxx = static_cast<float>(v[4 + np]);
            r.mzz = static_cast<float>(v[6 + np]);
            r.mxz = static_cast<float>(v[8 + np]);
            if (fmt[2] == 'm' && fmt[3] == 'w') {
                require(m <= 11.f, "magnitude must be <= 11");
                r.mo = seismic_moment(m);
            } else {
                r.mo = m;
            }
            if (ll) {
                check_ll(static_cast<float>(v[0]), static_cast<float>(v[1]));
                geomap_g2c(static_cast<float>(v[0]), static_cast<float>(v[1]), s.clon, s.clat, s.phi, r.sx, r.sy);
            } else {
                r.sx = static_cast<float>(v[0]);
                r.sy = static_cast<float>(v[1]);
            }
            r.sz = static_cast<float>(v[2]);
        } else if (fmt == "xym0dc" || fmt == "llm0dc" || fmt == "xymwdc" || fmt == "llmwdc") {
            auto v = read_reals(line, 3 + np + 4);
            take_prm(v, 3);
            const float m = static_cast<float>(v[3 + np]);
            const float strike = static_cast<float>(v[4 + np]), dip = static_cast<float>(v[5 + np]),
                        rake = static_cast<float>(v[6 + np]);
            check_sdr(strike, dip, rake);
            if (fmt[3] == 'w') {
                require(m <= 11.f, "magnitude must be <= 11");
                r.mo = seismic_moment(m);
            } else {
                r.mo = m;
            }
            sdr2moment(strike - s.phi, dip, rake, r.mxx, rd1, r.mzz, rd2, r.mxz, rd3);
            if (ll) {
                check_ll(static_cast<float>(v[0]), static_cast<float>(v[1]));
                geomap_g2c(static_cast<float>(v[0]), static_cast<float>(v[1]), s.clon, s.clat, s.phi, r.sx, r.sy);
            } else {
                r.sx = static_cast<float>(v[0]);
                r.sy = static_cast<float>(v[1]);
            }
            r.sz = static_cast<float>(v[2]);
        } else if (fmt == "xydsdc" || fmt == "lldsdc") {
            auto v = read_reals(line, 3 + np + 5);
            take_prm(v, 3);
            const float D = static_cast<float>(v[3 + np]), S = static_cast<float>(v[4 + np]);
            const float strike = static_cast<float>(v[5 + np]), dip = static_cast<float>(v[6 + np]),
                        rake = static_cast<float>(v[7 + np]);
            check_sdr(strike, dip, rake);
            sdr2moment(strike - s.phi, dip, rake, r.mxx, rd1, r.mzz, rd2, r.mxz, rd3);
            if (ll) {
                check_ll(static_cast<float>(v[0]), static_cast<float>(v[1]));
                geomap_g2c(static_cast<float>(v[0]), static_cast<float>(v[1]), s.clon, s.clat, s.phi, r.sx, r.sy);
            } else {
                r.sx = static_cast<float>(v[0]);
                r.sy = static_cast<float>(v[1]);
            }
            r.sz = static_cast<float>(v[2]);
            const int is0 = x2i(r.sx, s.xbeg, static_cast<float>(s.dx));
            const int ks0 = ef ? z2k(static_cast<float>(-R_EARTH * std::log((R_EARTH - r.sz) / R_EARTH)), s.zbeg,
                                     static_cast<float>(s.dz))
                               : z2k(r.sz, s.zbeg, static_cast<float>(s.dz));
            if (s.ibeg - 2 <= is0 && is0 <= s.iend + 3 && s.kbeg - 2 <= ks0 && ks0 <= s.kend + 3)
                r.mo = (1e9f * s.mu(ks0, is0)) * D * S;
            else
                r.mo = 0.f;
        } else if (fmt == "psmeca") {
            auto v = read_reals(line, 10);
            check_ll(static_cast<float>(v[0]), static_cast<float>(v[1]));
            geomap_g2c(static_cast<float>(v[0]), static_cast<float>(v[1]), s.clon, s.clat, s.phi, r.sx, r.sy);
            r.sz = static_cast<float>(v[2]);
            r.mzz = static_cast<float>(v[3]);
            r.mxx = static_cast<float>(v[4]);
            r.mxz = static_cast<float>(v[6]);
            const int iex = static_cast<int>(v[9]);
            const double M0tmp = std::sqrt(r.mxx * r.mxx + r.mzz * r.mzz + 2 * (r.mxz * r.mxz)) / std::sqrt(2.0f);
            r.mo = static_cast<float>(M0tmp * powi_sp(10.f, iex));
            r.prm[0] = 0.0f;
            const float half_dur_coef = 2 * 1.05f * 1e-8f;  // single-precision constant, as in Fortran
            r.prm[1] = static_cast<float>(half_dur_coef * std::pow(static_cast<double>(r.mo), 1.0 / 3.0));
            r.mo = r.mo * 1e-7f;
            r.mxx = static_cast<float>(r.mxx / M0tmp);
            r.mzz = static_cast<float>(r.mzz / M0tmp);
            r.mxz = static_cast<float>(r.mxz / M0tmp);
        } else {
            throw std::runtime_error("source: invalid stf_format '" + fmt + "'");
        }
        if (n == 0) {
            geomap_c2g(r.sx, r.sy, s.clon, s.clat, s.phi, s.evlo, s.evla);
            s.sx0 = r.sx; s.sy0 = 0; s.evdp = r.sz;
            s.mxx0 = r.mxx; s.myy0 = -12345.0f; s.mzz0 = r.mzz;
            s.myz0 = -12345.0f; s.mxz0 = r.mxz; s.mxy0 = -12345.0f;
            s.otim = r.prm[0];
        }
        out.push_back(r);
    }
}

void grid_bodyforce(State& s, const std::string& fn, const std::string& fmt, int np, std::vector<RawSource>& out) {
    const auto lines = data_lines(fn);
    const std::string coord = fmt.substr(0, 2);
    for (size_t n = 0; n < lines.size(); ++n) {
        RawSource r;
        r.prm.assign(np, 0.0f);
        auto v = read_reals(lines[n], 3 + np + 3);
        for (int j = 0; j < np; ++j) r.prm[j] = static_cast<float>(v[3 + j]);
        r.fx = static_cast<float>(v[3 + np]);
        r.fz = static_cast<float>(v[5 + np]);
        r.sz = static_cast<float>(v[2]);
        if (coord == "xy") {
            r.sx = static_cast<float>(v[0]);
        } else if (coord == "ll") {
            check_ll(static_cast<float>(v[0]), static_cast<float>(v[1]));
            float dum;
            geomap_g2c(static_cast<float>(v[0]), static_cast<float>(v[1]), s.clon, s.clat, s.phi, r.sx, dum);
        } else {
            throw std::runtime_error("source: invalid stf_format '" + fmt + "'");
        }
        if (n == 0) {
            geomap_c2g(r.sx, 0.0f, s.clon, s.clat, s.phi, s.evlo, s.evla);
            s.sx0 = r.sx; s.sy0 = 0.0f; s.evdp = r.sz; s.otim = r.prm[0];
            s.fx0 = r.fx; s.fy0 = -12345.0f; s.fz0 = r.fz;
        }
        out.push_back(r);
    }
}
}  // namespace

void source_setup(State& s, Source& src, const IniFile& prm) {
    s.pw_mode = prm.get_bool("pw_mode", false);
    if (s.pw_mode && !s.benchmark_mode)
        throw std::runtime_error("pw_mode (plane-wave incidence) is outside the scope of the C++ port");
    s.bf_mode = prm.get_bool("bf_mode", false);
    const std::string fn_stf = prm.get_string("fn_stf", "");
    src.stftype = prm.get_string("stftype", "kupper");
    const std::string stf_format = prm.get_string("stf_format", "xym0ij");
    src.sdep_fit = prm.get_string("sdep_fit", "asis");
    const bool ef = prm.get_bool("earth_flattening", false);
    if (src.stftype == "scosine") src.stftype = "cosine";
    src.n_stfprm = (src.stftype == "asymcos") ? 3 : 2;

    std::vector<RawSource> raw;
    if (s.benchmark_mode) {
        src.stftype = "kupper";
        src.n_stfprm = 2;
        s.bf_mode = false;
        RawSource r;
        r.sx = 0.0f; r.sz = 5.0f; r.mo = 1e15f;
        r.mxx = 1 / std::sqrt(2.0f); r.mzz = 1 / std::sqrt(2.0f); r.mxz = 0.0f;
        r.prm = {0.1f, 2.0f};
        raw.push_back(r);
    } else if (s.bf_mode) {
        grid_bodyforce(s, fn_stf, stf_format, src.n_stfprm, raw);
    } else {
        grid_moment(s, fn_stf, stf_format, src.n_stfprm, ef, raw);
    }
    require(!raw.empty(), "no sources defined");

    if (ef)
        for (auto& r : raw) r.sz = static_cast<float>(-R_EARTH * std::log((R_EARTH - r.sz) / R_EARTH));

    if (s.bf_mode) {
        float ss = 0.0f;
        for (auto& r : raw) ss += r.fx * r.fx + r.fz * r.fz;
        s.M0 = std::sqrt(ss);
        s.UC = s.UC * 1000.0f;
    } else {
        float ss = 0.0f;
        for (auto& r : raw) ss += r.mo;
        s.M0 = ss;
    }
    s.fcut = 0.f;
    for (auto& r : raw) s.fcut = std::max(s.fcut, 1 / r.prm[1]);
    s.fmax = 2 * s.fcut;

    for (auto& r : raw) {
        const int is = x2i(r.sx, s.xbeg, static_cast<float>(s.dx));
        int ks = z2k(r.sz, s.zbeg, static_cast<float>(s.dz));
        if (!(s.ibeg - 2 <= is && is <= s.iend + 3 && s.kbeg - 2 <= ks && ks <= s.kend + 3)) continue;
        float sz = r.sz;
        for (int k = 0; k <= NBD; ++k) {
            if (src.sdep_fit == "bd" + std::to_string(k)) {
                sz = s.bddep(k, is);
                ks = z2k(sz, s.zbeg, static_cast<float>(s.dz));
            }
        }
        require(s.xbeg <= r.sx && r.sx <= s.xend, "source x outside the model");
        require(s.zbeg <= sz && sz <= s.zend, "source z outside the model");
        src.isrc.push_back(is);
        src.ksrc.push_back(ks);
        src.sx.push_back(r.sx);
        src.sz.push_back(sz);
        src.srcprm.push_back(r.prm);
        if (s.bf_mode) {
            src.fx.push_back(static_cast<real_mp>(r.fx) / s.M0);
            src.fz.push_back(static_cast<real_mp>(r.fz) / s.M0);
        } else {
            src.mo.push_back(static_cast<real_mp>(r.mo) / s.M0);
            src.mxx.push_back(r.mxx);
            src.mzz.push_back(r.mzz);
            src.mxz.push_back(r.mxz);
        }
    }
    src.nsrc = static_cast<int>(src.isrc.size());
    src.dt_dxz = s.dt / (s.dx * s.dz);
}

void source_stressglut(State& s, const Source& src, int it) {
    if (s.bf_mode) return;
    for (int i = 0; i < src.nsrc; ++i) {
        const real_mp stime = momentrate(s.tbeg + (it - 0.5f) * s.dt, src.stftype, src.srcprm[i].data());
        const real_mp sdrop = src.mo[i] * stime * src.dt_dxz;
        const int ii = src.isrc[i], kk = src.ksrc[i];
        s.Sxx(kk, ii) -= src.mxx[i] * sdrop;
        s.Szz(kk, ii) -= src.mzz[i] * sdrop;
        s.Sxz(kk, ii) -= src.mxz[i] * sdrop / 4;
        s.Sxz(kk - 1, ii) -= src.mxz[i] * sdrop / 4;
        s.Sxz(kk, ii - 1) -= src.mxz[i] * sdrop / 4;
        s.Sxz(kk - 1, ii - 1) -= src.mxz[i] * sdrop / 4;
    }
}

void source_bodyforce(State& s, const Source& src, int it) {
    if (!s.bf_mode) return;
    for (int i = 0; i < src.nsrc; ++i) {
        const real_mp stime = momentrate(s.tbeg + it * s.dt, src.stftype, src.srcprm[i].data());
        const int ii = src.isrc[i], kk = src.ksrc[i];
        const real_mp r = s.rho(kk, ii);
        s.Vx(kk, ii) = static_cast<float>(s.Vx(kk, ii) + src.fx[i] / r * stime * src.dt_dxz / 2);
        s.Vx(kk, ii - 1) = static_cast<float>(s.Vx(kk, ii - 1) + src.fx[i] / r * stime * src.dt_dxz / 2);
        s.Vz(kk, ii) = static_cast<float>(s.Vz(kk, ii) + src.fz[i] / r * stime * src.dt_dxz / 2);
        s.Vz(kk - 1, ii) = static_cast<float>(s.Vz(kk - 1, ii) + src.fz[i] / r * stime * src.dt_dxz / 2);
    }
}

}  // namespace swpc
