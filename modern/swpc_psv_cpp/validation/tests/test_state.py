"""Full-wavefield (state-level) comparison: every grid cell of Vx, Vz, Sxx, Szz, Sxz and the
medium arrays must be bit-identical between the instrumented Fortran reference and the port."""
import os

import numpy as np
import pytest

from swpcval import runner, statedump as sd
from swpcval.scenario import Scenario

CASES = [
    Scenario("state_pml_visco", qp0=80.0, qs0=40.0, nt=300),
    Scenario("state_cerjan", abc_type="cerjan", nt=300),
    Scenario("state_lhm_ocean", vmodel_type="lhm", is_ocean=True, munk_profile=True, zbeg=-1.5, nt=300,
             lhm=" -1.0 1.0 1.5 0.0 1e6 1e6\n 0.3 2.2 3.8 2.0 150 80\n 2.0 2.8 6.0 3.5 600 300\n"),
]


@pytest.mark.differential
@pytest.mark.parametrize("sc", CASES, ids=[c.name for c in CASES])
def test_wavefield_bitwise_every_25_steps(sc, workdir, cpp_bin, record, monkeypatch):
    if not runner.fortran_bin("instrumented").exists():
        pytest.fail("instrumented Fortran build missing (build_fortran_ref.sh)")
    dumps = {}
    for eng in ("fortran", "cpp"):
        d = workdir / f"dump_{eng}"
        d.mkdir(parents=True, exist_ok=True)
        for f in d.glob("*.bin"):
            f.unlink()
        monkeypatch.setenv("SWPC_DUMP_DIR", str(d))
        monkeypatch.setenv("SWPC_DUMP_EVERY", "25")
        runner.run(eng, sc, workdir / eng, flavor="instrumented")
        dumps[eng] = d
    mf, mc = sd.read_medium(dumps["fortran"] / "medium.bin"), sd.read_medium(dumps["cpp"] / "medium.bin")
    assert mf["bounds"] == mc["bounds"]
    for f in sd.MEDIUM_FIELDS:
        assert np.array_equal(mf[f], mc[f]), f"medium {f} differs"
    steps = sd.steps(dumps["fortran"])
    assert steps == sd.steps(dumps["cpp"]) and len(steps) == sc.nt // 25
    cells = 0
    for it in steps:
        a = sd.read_state(dumps["fortran"] / f"state_{it:08d}.bin")
        b = sd.read_state(dumps["cpp"] / f"state_{it:08d}.bin")
        for f, _ in sd.STATE_FIELDS:
            n = int(np.count_nonzero(a[f] != b[f]))
            assert n == 0, f"step {it}: {f} differs in {n} cells"
            cells += a[f].size
        assert np.max(np.abs(a["Vx"])) > 0
    record(scenario=sc.name, checkpoints=len(steps), cell_values_compared=cells, bitwise=True)
