"""The OpenMP-parallel port must be bit-for-bit independent of the thread count."""
import pytest

from swpcval import runner
from swpcval.scenario import Scenario


@pytest.mark.parametrize("abc", ["pml", "cerjan"])
def test_thread_count_independence(abc, workdir, cpp_bin, record):
    sc = Scenario(f"det_{abc}", abc_type=abc, vmodel_type="lhm", nx=240, nz=140, nt=500)
    ref = None
    for n in (1, 2, 4, 8):
        runner.run("cpp", sc, workdir / f"t{n}", threads=n)
        out = {k: v.data.tobytes() for k, v in runner.outputs(workdir / f"t{n}").items()}
        if ref is None:
            ref = out
        assert out == ref, f"{n} threads differ from 1 thread"
    record(abc=abc, threads=[1, 2, 4, 8], files=len(ref), bitwise=True)


def test_repeat_runs_identical(workdir, cpp_bin):
    sc = Scenario("det_repeat", nt=200)
    outs = []
    for r in range(3):
        runner.run("cpp", sc, workdir / f"r{r}", threads=4)
        outs.append({k: v.data.tobytes() for k, v in runner.outputs(workdir / f"r{r}").items()})
    assert outs[0] == outs[1] == outs[2]
