// Port of the SAC writer in src/shared/m_sac.f90 (binary SAC, little-endian native).
#pragma once
#include <string>
#include <vector>

namespace swpc {

struct SacHeader {
    float f[70];          // words 1..70 (index = word - 1)
    int i[35];            // words 71..105
    int l[5];             // words 106..110
    std::string k[23];    // kstnm, kevnm(16), khole, ko, ka, kt0..kt9, kf, kuser0..2, kcmpnm, knetwk, kdatrd, kinst
    // npts, depmin/depmax/depmen and e are written exactly as set (the Fortran writer does not derive them).
    SacHeader();          // == sac__init
    float& delta() { return f[0]; }
    // named setters for the fields used by swpc_psv
};

enum SacF { DELTA = 1, DEPMIN = 2, DEPMAX = 3, B = 6, E = 7, O = 8, STLA = 32, STLO = 33, STDP = 35, EVLA = 36,
            EVLO = 37, EVDP = 39, MAG = 40, USER0 = 41, USER2 = 43, USER4 = 45, USER6 = 47, USER7 = 48,
            USER8 = 49, DIST = 51, AZ = 52, BAZ = 53, DEPMEN = 57, CMPAZ = 58, CMPINC = 59 };
enum SacI { NZYEAR = 71, NZJDAY = 72, NZHOUR = 73, NZMIN = 74, NZSEC = 75, NZMSEC = 76, NVHDR = 77, NPTS = 80,
            IFTYPE = 86, IDEP = 87 };
enum SacL { LEVEN = 106, LPSPOL = 107, LOVROK = 108, LCALDA = 109 };
enum SacK { KSTNM = 0, KEVNM = 1, KCMPNM = 19 };

void sac_write(const std::string& fn, const SacHeader& sh, const std::vector<float>& dat);

}  // namespace swpc
