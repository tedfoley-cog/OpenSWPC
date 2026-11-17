// Port of src/shared/m_readini.f90: "key = value" parameter files.
// First matching line wins, '#'/'!' lines are comments, values are parsed with
// Fortran list-directed rules (first token, quotes optional), ${VAR} expanded.
#pragma once
#include <string>
#include <vector>

namespace swpc {

class IniFile {
public:
    explicit IniFile(const std::string& path);
    bool has(const std::string& key) const;
    std::string get_string(const std::string& key, const std::string& def) const;
    double get_double(const std::string& key, double def) const;
    float get_float(const std::string& key, float def) const;
    int get_int(const std::string& key, int def) const;
    bool get_bool(const std::string& key, bool def) const;
    const std::string& path() const { return path_; }

private:
    bool find_raw(const std::string& key, std::string& rest) const;
    std::string path_;
    std::vector<std::string> lines_;
};

std::string list_directed_token(const std::string& s);
std::string expand_env(const std::string& s);
double parse_fortran_real(const std::string& tok);
bool parse_fortran_logical(const std::string& tok);

}  // namespace swpc
