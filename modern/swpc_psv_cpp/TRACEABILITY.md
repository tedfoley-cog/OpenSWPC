# Traceability: Fortran routine → C++ function → validation

Fortran paths are relative to `src/`; C++ paths to `modern/swpc_psv_cpp/src/`.
Test ids: **D** `test_differential` (end-to-end SAC, bitwise), **S** `test_state` (full wavefield,
bitwise), **R** `test_routines` (routine probes, bitwise), **T** `test_determinism`,
**E** `test_edge_cases`, **P** `test_physics`.

| Fortran | C++ | Validated by |
|---|---|---|
| `swpc_psv/main.f90` time loop | `main.cpp` `main` | D, S, T |
| `m_global` `global__readprm`, `global__setup`, `global__setup2`, `set_box`, `set_mpi_table` | `global.cpp` `global_setup`, `global_setup2`, `set_box` (single rank) | D, S (medium + bounds) |
| `m_global` `global__comm_vel`, `global__comm_stress`, `inside_node` | not needed (single process) | — |
| `m_kernel` `kernel__setup`, `kernel__update_stress`, `kernel__update_vel`, `kernel__vmax` | `kernel.cpp` `kernel_setup`, `kernel_update_stress`, `kernel_update_vel`, `make_fd_coef`, `order_sign`, `harmonic_mu` | S, D, T, P (analytic, convergence) |
| `m_absorb` dispatch | `kernel.cpp` `absorb_setup`, `absorb_update_stress`, `absorb_update_vel` | D, S |
| `m_absorb_p` `absorb_p__setup`, `damping_profile`, `update_vel_core`, `update_stress_core`, `set_*_boundary` | `kernel.cpp` `pml_setup`, `damping_profile`, `pml_update_vel_core`, `pml_update_stress_core` | R (`pml`), S, D, P (boundary reflection) |
| `m_absorb_c` `absorb_c__setup`, `absorb_c__update_*` | `kernel.cpp` `cerjan_setup`, `absorb_update_*` | S (`state_cerjan`), D, P |
| `m_medium` `medium__setup`, `relaxed_medium`, `surface_detection`, `velocity_minmax`, `stabilize_absorber` | `medium.cpp` `medium_setup`, `set_cell`, `for_all_media`, `stabilize_absorber` | S (medium arrays), D (`lhm_*`) |
| `m_vmodel_uni` `vmodel_uni` | `medium.cpp` `vmodel_uni` | D, S |
| `m_vmodel_lhm` `vmodel_lhm` | `medium.cpp` `vmodel_lhm`, `flattening` | D (`lhm_*`), S (`state_lhm_ocean`) |
| `shared/m_seawater` `seawater__init`, `seawater__vel` | `medium.cpp` `seawater_vel` | R (`munk`), D (`lhm_ocean_*`) |
| `m_source` `source__setup`, `source__grid_moment`, `source__grid_bodyforce` | `source.cpp` `source_setup`, `grid_moment`, `grid_bodyforce`, `check_sdr`, `check_ll` | D (`src_*`), E |
| `m_source` `source__stressglut`, `source__bodyforce` | `source.cpp` `source_stressglut`, `source_bodyforce` | D, S, P (reciprocity) |
| `m_source` `pw_setup` | not migrated (rejected) | E |
| `shared/m_fdtool` source time functions, `momentrate` | `fdtool.cpp` `kupper`, `texp`, `cosine`, `boxcar`, `triangle`, `herrmann`, `asymcos`, `momentrate` | R (`stf`, 2,807 cases), D |
| `shared/m_fdtool` `sdr2moment`, `seismic_moment`, `moment_magnitude` | `fdtool.cpp` same names | R (`sdr`, `mw`) |
| `shared/m_fdtool` `x2i`, `z2k`, `i2x`, `k2z` | `fdtool.cpp` same names | R (`grid`) |
| `shared/m_fdtool` `visco_set_relaxtime`, `visco_chi`, `visco_constq_zeta` | `fdtool.cpp` same names | R (`visco`), P (Q) |
| `shared/m_fdtool` stability helpers; `m_report` `report__setup` CFL check | `fdtool.cpp` `fdm_cond_stability`, `fdm_cond_wavelength`; `main.cpp` `report_setup` | E (CFL reject / accept) |
| `shared/m_geomap`, `shared/m_gk` | `geomap.cpp` `geomap_g2c`, `geomap_c2g`, `lltoxy`, `xytoll` | R (`geo`) |
| `shared/m_readini` | `readini.cpp` `IniFile`, `parse_fortran_real`, `parse_fortran_logical` | E (parser semantics), D (`input_fortran_syntax`) |
| `m_wav` `wav__setup`, `set_stinfo`, `wav__store` | `wav.cpp` `wav_setup`, `set_stinfo`, `wav_store` | D (all components) |
| `m_wav` `wav__write`, `set_sac_header`, `initialize_sac_header`, `export_wav__sac`; `shared/m_sac` | `wav_write.cpp` `wav_write`, `base_header`; `sac.cpp` `sac_write` | D (header word-by-word + samples) |
| `m_wav` `export_wav__csf`, `export_wav__tar`; `m_snap` (all) | not migrated (rejected) | E |
| `m_vmodel_grd*`, `m_vmodel_lgm*`, `m_vmodel_*_rmed`, `m_vmodel_user` | not migrated (rejected) | E |
| — | `dump.cpp` `dump_medium`, `dump_state` (validation instrumentation; Fortran twin in `validation/fortran_instrument/dump.f90.inc`) | S |
