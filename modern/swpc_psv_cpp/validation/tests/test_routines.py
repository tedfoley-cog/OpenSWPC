"""Routine-level differential test: source time functions, moment tensors, magnitude, grid index
maps, viscoelastic relaxation, map projection, seawater profile and PML profiles are evaluated
by upstream Fortran objects and by the C++ port on identical inputs and compared bit-for-bit."""
import collections
import os
import subprocess

import pytest

from swpcval import runner

HERE = runner.ROOT / "validation" / "routines"
OUT = runner.ROOT / "build" / "probe"


@pytest.mark.differential
def test_routines_bitwise(record):
    subprocess.run([str(HERE / "build_probe.sh")], check=True, capture_output=True)
    env = dict(os.environ, HWLOC_COMPONENTS="no_os,stop")
    f = subprocess.run([str(OUT / "probe_fortran")], check=True, capture_output=True, text=True, env=env).stdout
    c = subprocess.run([str(OUT / "probe_cpp")], check=True, capture_output=True, text=True).stdout
    fl = [" ".join(x.split()) for x in f.splitlines() if x.strip()]
    cl = [" ".join(x.split()) for x in c.splitlines() if x.strip()]
    assert len(fl) == len(cl) > 8000
    tot, bad = collections.Counter(), collections.Counter()
    first = {}
    for a, b in zip(fl, cl):
        tag = a.split()[0]
        tot[tag] += 1
        if a != b:
            bad[tag] += 1
            first.setdefault(tag, (a, b))
    record(cases=dict(tot), mismatches=dict(bad))
    assert not bad, f"routine mismatches: {dict(bad)}; first: {first}"
