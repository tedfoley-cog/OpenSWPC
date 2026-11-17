#!/usr/bin/env bash
# Build the upstream Fortran swpc_psv.x as the differential-test reference.
#   strict  : -O2 without -ffast-math and with FMA contraction disabled (IEEE-reproducible reference)
#   instrumented: strict + full-wavefield dumps (SWPC_DUMP_DIR / SWPC_DUMP_EVERY) for state-level comparison
#   upstream: the stock ubuntu-gfortran flags (-O2 -ffast-math), to measure what the shipped binary does
# Usage: build_fortran_ref.sh [outdir]   (default: ../build/fortran-ref)
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
repo="$(cd "$here/../../.." && pwd)"
out="${1:-$here/../build/fortran-ref}"
mkdir -p "$out"
for flavor in strict upstream instrumented; do
  work="$out/src-$flavor"
  rm -rf "$work"
  mkdir -p "$work/bin"
  cp -r "$repo/src" "$work/src"
  if [ "$flavor" = instrumented ]; then
    m="$work/src/swpc_psv/main.f90"
    python3 - "$m" "$here/fortran_instrument/dump.f90.inc" <<'PY'
import sys
m, inc = sys.argv[1], sys.argv[2]
s = open(m).read()
s = s.replace("    call absorb__setup(io_prm)\n", "    call absorb__setup(io_prm)\n    call vdump_medium()\n", 1)
s = s.replace("        call global__comm_vel()\n", "        call global__comm_vel()\n        call vdump_state(it)\n", 1)
s = s.replace("end program swpc_psv", open(inc).read() + "end program swpc_psv", 1)
assert s.count("vdump_") >= 6, "instrumentation patch did not apply"
open(m, "w").write(s)
PY
  fi
  if [ "$flavor" != upstream ]; then
    flags="-O2 -fopenmp -cpp -ffp-contract=off"
  else
    flags="-O2 -ffast-math -fopenmp -cpp"
  fi
  make -C "$work/src/swpc_psv" arch=ubuntu-gfortran FFLAGS="$flags" -j"$(nproc)" >"$out/build-$flavor.log" 2>&1
  cp "$work/bin/swpc_psv.x" "$out/swpc_psv_$flavor.x"
  echo "built $out/swpc_psv_$flavor.x ($flags)"
done
