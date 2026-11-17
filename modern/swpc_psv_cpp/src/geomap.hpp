// Port of m_geomap.f90 / m_gk.f90 (Gauss-Krueger projection around clon/clat).
#pragma once
namespace swpc {
void geomap_g2c(float lon, float lat, float lon0, float lat0, float phi, float& x, float& y);
void geomap_c2g(float x, float y, float lon0, float lat0, float phi, float& lon, float& lat);
}  // namespace swpc
