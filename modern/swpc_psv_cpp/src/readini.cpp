#include "readini.hpp"

#include <cstdlib>
#include <fstream>
#include <stdexcept>

namespace swpc {

namespace {
std::string adjustl(const std::string& s) {
    size_t p = s.find_first_not_of(" \t");
    return p == std::string::npos ? std::string() : s.substr(p);
}
std::string trim(const std::string& s) {
    std::string a = adjustl(s);
    size_t p = a.find_last_not_of(" \t\r\n");
    return p == std::string::npos ? std::string() : a.substr(0, p + 1);
}
}  // namespace

IniFile::IniFile(const std::string& path) : path_(path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open parameter file: " + path);
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines_.push_back(line);
    }
}

bool IniFile::find_raw(const std::string& key, std::string& rest) const {
    const std::string keyword = trim(key);
    for (const auto& raw : lines_) {
        std::string cline = adjustl(raw);
        if (cline.empty() || cline[0] == '#' || cline[0] == '!') continue;
        if (cline.compare(0, keyword.size(), keyword) != 0) continue;
        std::string r = adjustl(cline.substr(keyword.size()));
        if (r.empty() || r[0] != '=') continue;
        rest = adjustl(r.substr(1));
        return true;
    }
    return false;
}

bool IniFile::has(const std::string& key) const {
    std::string r;
    return find_raw(key, r);
}

std::string list_directed_token(const std::string& s0) {
    std::string s = adjustl(s0);
    if (s.empty()) return "";
    if (s[0] == '\'' || s[0] == '"') {
        const char q = s[0];
        std::string out;
        for (size_t i = 1; i < s.size(); ++i) {
            if (s[i] == q) {
                if (i + 1 < s.size() && s[i + 1] == q) {
                    out += q;
                    ++i;
                } else {
                    break;
                }
            } else {
                out += s[i];
            }
        }
        return out;
    }
    size_t e = s.find_first_of(" \t,/");
    return e == std::string::npos ? s : s.substr(0, e);
}

std::string expand_env(const std::string& s) {
    std::string out;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '$' && i + 1 < s.size() && s[i + 1] == '{') {
            size_t e = s.find('}', i + 2);
            if (e != std::string::npos) {
                const char* v = std::getenv(s.substr(i + 2, e - i - 2).c_str());
                if (v) out += v;
                i = e;
                continue;
            }
        }
        out += s[i];
    }
    return out;
}

double parse_fortran_real(const std::string& tok) {
    std::string t = tok;
    for (auto& c : t)
        if (c == 'd' || c == 'D' || c == 'q' || c == 'Q') c = 'e';
    size_t pos = 0;
    double v = std::stod(t, &pos);
    if (pos != t.size()) throw std::invalid_argument("bad real: " + tok);
    return v;
}

bool parse_fortran_logical(const std::string& tok) {
    std::string t = tok;
    if (!t.empty() && t[0] == '.') t = t.substr(1);
    if (t.empty()) throw std::invalid_argument("bad logical: " + tok);
    if (t[0] == 't' || t[0] == 'T') return true;
    if (t[0] == 'f' || t[0] == 'F') return false;
    throw std::invalid_argument("bad logical: " + tok);
}

std::string IniFile::get_string(const std::string& key, const std::string& def) const {
    std::string rest;
    if (!find_raw(key, rest)) return def;
    return expand_env(list_directed_token(rest));
}

double IniFile::get_double(const std::string& key, double def) const {
    std::string rest;
    if (!find_raw(key, rest)) return def;
    return parse_fortran_real(list_directed_token(rest));
}

float IniFile::get_float(const std::string& key, float def) const {
    std::string rest;
    if (!find_raw(key, rest)) return def;
    return std::stof([&] {
        std::string t = list_directed_token(rest);
        for (auto& c : t)
            if (c == 'd' || c == 'D') c = 'e';
        return t;
    }());
}

int IniFile::get_int(const std::string& key, int def) const {
    std::string rest;
    if (!find_raw(key, rest)) return def;
    const std::string t = list_directed_token(rest);
    size_t pos = 0;
    int v = std::stoi(t, &pos);
    if (pos != t.size()) throw std::invalid_argument("bad integer for " + key + ": " + t);
    return v;
}

bool IniFile::get_bool(const std::string& key, bool def) const {
    std::string rest;
    if (!find_raw(key, rest)) return def;
    return parse_fortran_logical(list_directed_token(rest));
}

}  // namespace swpc
