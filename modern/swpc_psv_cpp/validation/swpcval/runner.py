"""Run the Fortran reference and the C++ port on the same scenario and collect SAC output."""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from .sac import Sac, read_sac
from .scenario import Scenario

ROOT = Path(__file__).resolve().parents[2]          # modern/swpc_psv_cpp
CPP_BIN = Path(os.environ.get("SWPC_CPP_BIN", ROOT / "build" / "swpc_psv_cpp"))
FORTRAN_DIR = Path(os.environ.get("SWPC_FORTRAN_REF_DIR", ROOT / "build" / "fortran-ref"))


def fortran_bin(flavor: str = "strict") -> Path:
    return FORTRAN_DIR / f"swpc_psv_{flavor}.x"


def run(engine: str, sc: Scenario, workdir: Path, threads: int = 1, flavor: str = "strict",
        check: bool = True) -> subprocess.CompletedProcess:
    """engine: 'cpp' or 'fortran'. Runs inside workdir so relative paths in the .inf resolve."""
    if workdir.exists():
        shutil.rmtree(workdir)
    inf = sc.write(workdir)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    if engine == "cpp":
        cmd = [str(CPP_BIN), "-i", inf.name, "-q"]
    elif engine == "fortran":
        cmd = ["mpirun", "--oversubscribe", "--allow-run-as-root", "-np", "1", str(fortran_bin(flavor)), "-i", inf.name]
        # hwloc topology discovery crashes on some virtualised hosts; MPI binding is irrelevant for -np 1.
        env.setdefault("HWLOC_COMPONENTS", "no_os,stop")
        # Serve every allocation from fresh zero-filled pages: the upstream code reads a few
        # uninitialised heap words (finding F-2 in VALIDATION_REPORT.md), which otherwise makes the
        # reference itself non-reproducible at the 1e-40 level.
        env.setdefault("MALLOC_MMAP_THRESHOLD_", "0")
    else:
        raise ValueError(engine)
    p = subprocess.run(cmd, cwd=workdir, env=env, capture_output=True, text=True, timeout=1800)
    if check and p.returncode != 0:
        raise RuntimeError(f"{engine} failed ({p.returncode}) on {sc.name}:\n{p.stdout[-3000:]}\n{p.stderr[-3000:]}")
    return p


def outputs(workdir: Path) -> dict[str, Sac]:
    """Map '<station>.<component>' -> Sac for every SAC file produced in workdir/out/wav."""
    res = {}
    for f in sorted((workdir / "out" / "wav").glob("*.sac")):
        s = read_sac(f)
        res[f"{s.kstnm}.{s.kcmpnm}"] = s
    return res
