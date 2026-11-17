# Validation report — swpc_psv C++17 migration

**Result: the C++ port reproduces the upstream Fortran solver bit-for-bit** on every scenario,
every waveform component and every grid cell checked. Physics checks that do not use the
Fortran code also pass. All numbers below come from `python3 -m pytest` in `validation/`
(49 tests, about 95 s on 8 cores) and are re-generated in CI by `validation/report.py`.

## 1. Approach

| Layer | What is compared | Acceptance |
|---|---|---|
| Routine probes (`test_routines.py`) | 8,300 evaluations of source time functions, moment tensors, Mw, grid index maps, relaxation times, map projection, seawater profile, PML profiles. Upstream Fortran objects vs C++ on identical inputs | bit-identical |
| Full wavefield (`test_state.py`) | Every cell of `Vx, Vz, Sxx, Szz, Sxz` every 25 steps, plus all medium arrays, from an instrumented Fortran build (scratch copy; upstream source not modified) | bit-identical |
| End-to-end (`test_differential.py`) | 23 scenarios × all SAC outputs (velocity, displacement, stress, strain); every header word except the wall-clock date | bit-identical samples, identical headers |
| Determinism (`test_determinism.py`) | C++ output at 1, 2, 4 and 8 OpenMP threads; repeated runs | bit-identical |
| Edge cases (`test_edge_cases.py`) | CFL rejection and acceptance, parser semantics, out-of-scope inputs, missing files | same accept/reject as Fortran; port fails loudly outside scope |
| Physics (`test_physics.py`) | Analytic solution, grid convergence, absorbing-boundary reflection, reciprocity, attenuation | thresholds in §3 |

The Fortran reference is built three ways by `validation/build_fortran_ref.sh`:

- **strict:** `-O2`, no `-ffast-math`, `-ffp-contract=off`. This is the comparison target.
- **instrumented:** strict plus wavefield dumps.
- **upstream:** the shipped `-ffast-math` flags. This build shows how much the upstream binary itself moves with compiler flags.

## 2. Fortran ↔ C++ equivalence

- **End-to-end:** 900 SAC files over 23 scenarios. All 900 are bit-identical (0 ULP) with identical headers. The scenarios cover:
  - uniform, layered, ocean (uniform water and Munk profile), topography, full space, Earth flattening and PML stabilisation models
  - PML and Cerjan absorbers
  - moment-tensor, strike/dip/rake, Mw and body-force sources
  - all 7 source time functions
  - anisotropic grids, a thin PML, a long run and Fortran `D`-exponent input syntax
- **Full wavefield:** 3 scenarios (viscoelastic PML; Cerjan; layered model with Munk ocean) × 12 checkpoints, 4.67 M cell values. All bit-identical.
- **Routines:** 8,300 / 8,300 probes bit-identical.
- **For scale:** the upstream `-ffast-math` build differs from the strict build by 1.45 × 10⁻⁶ relative L∞. The C++ port sits at the same distance from that build because it is identical to the strict build.

## 3. Physics checks (independent of the Fortran code)

| Check | Setup | Result | Threshold |
|---|---|---|---|
| Analytic solution | Explosion line source in a full space; 5 receivers at 2–4.3 km; 2-D Green's function via \(H^{(2)}_1\); one amplitude scale fitted across all traces | 0.63 % L2 misfit | < 3 % |
| Grid convergence | dx = 100, 50, 25 m against the analytic solution | 2.76 % → 1.00 % → 0.63 %; order 1.47 on the coarse pair, then flattening at the comparison floor | monotone; coarse order > 1.2 |
| PML reflection | 10 km domain vs 32 km reference domain, 20-cell PML | 0.055 % of peak | < 2 % |
| Cerjan reflection | same, 20-cell sponge | 5.8 % of peak | < 15 % |
| Reciprocity | Layered, attenuating model; \(G_{xz}(B;A)\) vs \(G_{zx}(A;B)\) | 0.90 % L2 | < 1 % |
| Constant Q | Q = 50 vs elastic run, spectral ratio between 2 and 5 km, 0.5–4 Hz | Q = 59.7 (53–62 depending on band) | ± 25 % |

The convergence and Q results are characteristics of the upstream numerical scheme, not of the port: the port is bit-identical to the Fortran code. Below 50 m, the misfit is set by how the point source is represented and how it is compared with the analytic solution, not by the grid. The Q estimate reflects the solver's approximation of constant Q with three relaxation mechanisms.

## 4. Findings

| ID | Type | Finding | Handling |
|---|---|---|---|
| F-1 | Upstream defect | `source__grid_bodyforce` never assigns `sy(i)` for `xy` input. The event lon/lat SAC header words (`evla`, `evlo`) are computed from uninitialised memory; NaN was observed | Port uses `sy = 0`. The differential test excludes those 2 header words for body-force runs only, and `test_known_defect_bodyforce_event_location` pins the correct value |
| F-2 | Upstream defect | The Fortran binary reads uninitialised heap memory. Values of about 10⁻⁴⁴ appeared at the domain edge from step 2 in Cerjan runs, so the reference is not reproducible run to run. The specific array was not isolated | Reference runs use `MALLOC_MMAP_THRESHOLD_=0` (fresh zero pages). With it, Cerjan runs are bit-identical; the port zero-initialises all arrays |
| F-3 | Port defect, caught and fixed | PML `damping_profile` evaluated `d0*r*r` instead of `d0*(r*r)`, a 1-ULP coefficient difference that grew to 1.8 × 10⁻⁵ in a topography run. Localised by wavefield bisection (first difference at step 66, PML cell) and confirmed by the routine probe | Fixed; `pml` probe cases guard it |
| F-4 | Port defect, caught and fixed | Cerjan profile had the same operation-order issue (`-alpha*r*r`) | Fixed |
| F-5 | Port defect, caught and fixed | SAC writer emitted 24 character-header fields instead of 23, giving the wrong file size and shifting `kcmpnm` | Fixed; header compared word by word |
| F-6 | Environment | OpenMPI's hwloc topology probe crashes on some VMs | Runner sets `HWLOC_COMPONENTS=no_os,stop` |

## 5. Scope and limits

The port covers the `uni` and `lhm` velocity models, the ocean column, PML and Cerjan absorbers, all source formats and source time functions, and SAC output, all in a single process (see `README.md`).

Out-of-scope inputs (`grd`, `lgm`, random media, user models, plane-wave mode, snapshots, `csf` and tar output) stop with an explicit error. This is tested.

MPI decomposition is not migrated. Equivalence is therefore shown against single-rank Fortran runs.
