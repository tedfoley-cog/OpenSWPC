# swpc_psv_cpp — C++17 migration of the OpenSWPC 2-D P-SV solver

A C++17 / OpenMP port of `src/swpc_psv` (the 2-D P-SV viscoelastic finite-difference
seismic-wave simulator) and the shared modules it depends on. It reads the same `.inf`
parameter files, source/station files and layered-model files as the Fortran code and writes the
same SAC waveforms.

The Fortran tree is untouched: it is the reference implementation that the port is validated
against. See [VALIDATION_REPORT.md](VALIDATION_REPORT.md) for the evidence and
[TRACEABILITY.md](TRACEABILITY.md) for the routine-by-routine map.

## Build and run

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
OMP_NUM_THREADS=8 build/swpc_psv_cpp -i path/to/input.inf
```

`-DSWPC_FAST_MATH=ON` adds `-ffast-math` (same as the upstream default Makefile). The default build
uses `-ffp-contract=off` so that results are bit-reproducible and comparable to the strict
Fortran reference.

## Scope of this migration

| Supported | Not yet migrated (the port stops with an error message) |
|---|---|
| `vmodel_type = uni, lhm`; ocean column (`is_ocean`, `munk_profile`); `earth_flattening`; `stabilize_pml` | `grd`, `lgm`, `user`, random-media variants |
| PML and Cerjan absorbers; `fullspace_mode`; topography from layered models | MPI domain decomposition (`nproc_x` accepted, run is single-process) |
| Moment-tensor, strike/dip/rake, Mw and body-force sources in `xy`/`ll` formats; all 7 source time functions | Plane-wave mode (`pw_mode`) |
| Stations in `xy`/`ll` with `dep`, `fsb`, `obb`, `oba`, `bd0`–`bd9`; velocity, displacement, stress, strain SAC output | NetCDF/`.snp` snapshots, `csf` waveforms, tar packing |

## Validation

```bash
validation/build_fortran_ref.sh          # strict, upstream-flag and instrumented Fortran builds
cd validation && python3 -m pytest -v    # ~2 minutes on 8 cores
python3 report.py                        # Markdown evidence tables from build/validation-results
```

Requires gfortran, OpenMPI, NetCDF-Fortran (to build the reference), CMake, and Python with
numpy, scipy and pytest. CI: `.github/workflows/swpc-psv-cpp-validation.yml`.

On hosts where OpenMPI's hwloc topology probe crashes, the runner already sets
`HWLOC_COMPONENTS=no_os,stop`.
