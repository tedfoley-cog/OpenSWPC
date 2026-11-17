#!/usr/bin/env bash
# Build the Fortran and C++ routine probes. Requires build_fortran_ref.sh to have run (strict flavor).
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"
ref="$root/build/fortran-ref/src-strict/src"
out="$root/build/probe"
mkdir -p "$out"
sed -n '/^ *subroutine damping_profile/,/^ *end subroutine damping_profile/p' "$ref/swpc_psv/m_absorb_p.f90" >"$out/pml_profile.inc"
grep -q 'end subroutine damping_profile' "$out/pml_profile.inc"
flags="-O2 -fopenmp -cpp -ffp-contract=off"
o="$ref/swpc_psv"
mpif90 $flags -I"$o" -I"$out" -J"$out" -o "$out/probe_fortran" "$here/probe.f90" \
  "$o/m_std.o" "$o/m_fdtool.o" "$o/m_gk.o" "$o/m_geomap.o" "$o/m_seawater.o"
g++ -std=c++17 -O2 -ffp-contract=off -fopenmp -I"$root/src" -o "$out/probe_cpp" "$here/probe.cpp" \
  "$root/src/fdtool.cpp" "$root/src/geomap.cpp" "$root/src/kernel.cpp" "$root/src/medium.cpp" "$root/src/readini.cpp"
echo "built $out/probe_fortran $out/probe_cpp"
