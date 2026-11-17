import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from swpcval import runner  # noqa: E402

RESULTS = runner.ROOT / "build" / "validation-results"
WORK = runner.ROOT / "build" / "validation-work"


def pytest_configure(config):
    RESULTS.mkdir(parents=True, exist_ok=True)
    if not hasattr(config, "workerinput"):
        for f in RESULTS.glob("*.jsonl"):
            f.unlink()


@pytest.fixture(scope="session")
def cpp_bin():
    if not runner.CPP_BIN.exists():
        pytest.fail(f"C++ binary not built: {runner.CPP_BIN} (cmake -S . -B build && cmake --build build)")
    return runner.CPP_BIN


@pytest.fixture(scope="session")
def fortran_ref():
    if not runner.fortran_bin("strict").exists():
        pytest.fail(f"Fortran reference not built: run validation/build_fortran_ref.sh")
    return runner.fortran_bin("strict")


@pytest.fixture
def workdir(request):
    d = WORK / request.node.name.replace("/", "_").replace("[", "-").replace("]", "")
    d.mkdir(parents=True, exist_ok=True)
    return d


@pytest.fixture
def record(request):
    """Persist structured evidence for the validation report: record(**fields)."""
    name = request.node.nodeid.split("::", 1)[1]
    group = Path(request.node.fspath).stem

    def _rec(**fields):
        f = RESULTS / f"{group}.jsonl"
        with f.open("a") as fh:
            fh.write(json.dumps({"test": name, **fields}, default=float) + "\n")

    return _rec
