"""Fortran (strict-IEEE build of upstream swpc_psv.x) vs C++ port on a matrix of scenarios.

Every SAC file the Fortran code writes must also be written by the C++ port, with an identical
header (except the wall-clock execution date) and BIT-IDENTICAL waveform samples. The float
tolerances below are only used to grade the size of a deviation in the failure message.
"""
import numpy as np
import pytest

from swpcval import runner
from swpcval.compare import compare
from swpcval.sac import header_diff
from swpcval.scenario import Scenario

# Grading thresholds (relative to the trace peak). float32 round-off accumulated over
# O(10^3) steps is ~1e-6; 1e-5 leaves headroom while catching any real algorithmic difference
# (a single wrong coefficient or index shows up at 1e-2 .. 1e0).
REL_LINF = 1e-5
REL_L2 = 1e-5

LAYERED_OCEAN = """# depth rho vp vs qp qs
 -1.0  1.0 1.5 0.0  1e6 1e6
  0.3  2.2 3.8 2.0  150 80
  1.5  2.5 5.0 2.9  300 150
  3.5  2.8 6.0 3.5  600 300
"""

DIP_SLIP = ["0.0 0.0 3.0  0.0 0.6  1.0e15  30.0 45.0 90.0"]
MW_DC = ["1.0 0.0 2.5  0.1 0.4  5.5  10.0 60.0 -90.0"]
FORCE = ["0.0 0.0 2.0  0.0 0.4  1.0e15  0.0 1.0e15"]
FORCE_X = ["0.5 0.0 2.0  0.0 0.4  1.0e15  0.0 0.0"]
MULTI = ["-1.0 0.0 2.0  0.0 0.5  1.0e15  1.0 0.0 -1.0 0.0 0.0 0.0",
         " 1.0 0.0 3.0  0.2 0.5  2.0e15  0.0 0.0 0.0 0.0 1.0 0.0",
         " 0.0 0.0 4.0  0.4 0.5  0.5e15  0.7 0.0 0.7 0.0 0.1 0.0"]

# Finding F-1: upstream source__grid_bodyforce never assigns sy(i) for 'xy' input, so the
# event lon/lat header words (36 evla, 37 evlo) are computed from uninitialised memory in the
# Fortran build. The port uses sy = 0 (the documented meaning); test_known_defects pins that.
KNOWN_UPSTREAM_HEADER_DEFECTS = {36, 37}

SCENARIOS = [
    Scenario("uni_elastic_pml"),
    Scenario("uni_visco_pml", qp0=80.0, qs0=40.0),
    Scenario("uni_cerjan", abc_type="cerjan"),
    Scenario("uni_fullspace", fullspace_mode=True, zbeg=-6.0, sources=["0.0 0.0 0.0  0.0 0.5  1.0e15  1.0 0.0 1.0 0.0 0.3 0.0"]),
    Scenario("uni_topo", topo0=0.5, zbeg=-1.0),
    Scenario("lhm_visco_pml", vmodel_type="lhm"),
    Scenario("lhm_ocean_munk", vmodel_type="lhm", lhm=LAYERED_OCEAN, is_ocean=True, munk_profile=True, zbeg=-1.5,
             stations=["-5.0 0.0 0.0 OB01 obb", "-2.0 0.0 0.0 OA02 oba", "2.5 0.0 1.5 ST03 dep", "6.0 0.0 0.0 FS04 fsb"]),
    Scenario("lhm_ocean_uniform_water", vmodel_type="lhm", lhm=LAYERED_OCEAN, is_ocean=True, munk_profile=False, zbeg=-1.5,
             stations=["-5.0 0.0 0.0 OB01 obb", "3.0 0.0 0.0 OB02 obb"]),
    Scenario("lhm_stabilize_pml", vmodel_type="lhm", stabilize_pml=True),
    Scenario("lhm_earth_flattening", vmodel_type="lhm", earth_flattening=True),
    Scenario("src_double_couple", stf_format="xym0dc", sources=DIP_SLIP),
    Scenario("src_mw_double_couple_texp", stf_format="xymwdc", stftype="texp", sources=MW_DC),
    Scenario("src_multi_cosine", stftype="cosine", sources=MULTI),
    Scenario("src_triangle", stftype="triangle"),
    Scenario("src_boxcar", stftype="boxcar"),
    Scenario("src_herrmann", stftype="herrmann"),
    Scenario("src_asymcos", stftype="asymcos", sources=["0.0 0.0 3.0  0.0 0.2 0.4  1.0e15  1.0 0.0 1.0 0.0 0.0 0.0"]),
    Scenario("src_bodyforce_z", bf_mode=True, stf_format="xy", sources=FORCE),
    Scenario("src_bodyforce_x_cerjan", bf_mode=True, stf_format="xy", sources=FORCE_X, abc_type="cerjan"),
    Scenario("grid_anisotropic", dx=0.1, dz=0.05, nz=200, dt=0.006),
    Scenario("grid_thin_pml", na=8, nt=500),
    Scenario("grid_long_run", nt=2000, nx=160, nz=100, ntdec_w=10),
    Scenario("input_fortran_syntax", extra={"dt": "8.0d-3", "vp0": "5.0D0", "qp0": "1.D6"}),
]


@pytest.mark.differential
@pytest.mark.parametrize("sc", SCENARIOS, ids=[s.name for s in SCENARIOS])
def test_fortran_vs_cpp(sc, workdir, cpp_bin, fortran_ref, record):
    runner.run("fortran", sc, workdir / "fortran")
    runner.run("cpp", sc, workdir / "cpp")
    fo = runner.outputs(workdir / "fortran")
    co = runner.outputs(workdir / "cpp")
    assert fo, "Fortran produced no SAC output"
    assert sorted(fo) == sorted(co), f"file sets differ: {sorted(set(fo) ^ set(co))}"

    ignore = KNOWN_UPSTREAM_HEADER_DEFECTS if sc.bf_mode else set()
    worst_linf = worst_l2 = 0.0
    nbit = 0
    max_ulp = 0
    for key in sorted(fo):
        hd = [d for d in header_diff(fo[key], co[key]) if int(d.split()[1]) not in ignore]
        assert not hd, f"{key}: SAC header mismatch: {hd}"
        m = compare(fo[key].data, co[key].data)
        assert np.all(np.isfinite(co[key].data)), f"{key}: non-finite samples"
        assert m.bitwise, (f"{key}: samples not bit-identical (rel L-inf {m.rel_linf:.3e}, rel L2 {m.rel_l2:.3e}, "
                           f"max {m.max_ulp} ULP; float32 round-off level would be < {REL_LINF})")
        worst_linf = max(worst_linf, m.rel_linf)
        worst_l2 = max(worst_l2, m.rel_l2)
        max_ulp = max(max_ulp, m.max_ulp)
        nbit += m.bitwise
    record(scenario=sc.name, files=len(fo), bitwise_files=nbit, worst_rel_linf=worst_linf, worst_rel_l2=worst_l2,
           max_ulp=max_ulp, headers_identical=True)


@pytest.mark.differential
def test_upstream_fast_math_build_within_tolerance(workdir, cpp_bin, record):
    """The stock upstream build uses -ffast-math; show the port is as close to it as it is to itself."""
    if not runner.fortran_bin("upstream").exists():
        pytest.fail("upstream Fortran build missing")
    sc = Scenario("upstream_fastmath", vmodel_type="lhm")
    runner.run("fortran", sc, workdir / "strict")
    runner.run("fortran", sc, workdir / "upstream", flavor="upstream")
    runner.run("cpp", sc, workdir / "cpp")
    st, up, cp = (runner.outputs(workdir / x) for x in ("strict", "upstream", "cpp"))
    fast_vs_strict = max(compare(st[k].data, up[k].data).rel_linf for k in st)
    cpp_vs_fast = max(compare(up[k].data, cp[k].data).rel_linf for k in st)
    record(fortran_fastmath_vs_strict=fast_vs_strict, cpp_vs_fortran_fastmath=cpp_vs_fast)
    assert cpp_vs_fast <= max(10 * fast_vs_strict, REL_LINF)


@pytest.mark.differential
def test_known_defect_bodyforce_event_location(workdir, cpp_bin, record):
    """F-1: for body-force sources the port writes the well-defined epicentre (c2g of (sx, 0))."""
    import numpy as np  # noqa: F811
    sc = Scenario("bf_evloc", bf_mode=True, stf_format="xy", sources=FORCE, nt=10)
    runner.run("cpp", sc, workdir / "cpp")
    co = runner.outputs(workdir / "cpp")
    h = next(iter(co.values()))
    # source at x = 0 on the map origin -> epicentre is (clon, clat)
    assert abs(h.f("evlo") - 139.7604) < 1e-4 and abs(h.f("evla") - 35.7182) < 1e-4
    assert np.isfinite(h.f("evlo")) and np.isfinite(h.f("evla"))
    record(evlo=h.f("evlo"), evla=h.f("evla"))
