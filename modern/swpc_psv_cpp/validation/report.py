#!/usr/bin/env python3
"""Render build/validation-results/*.jsonl (written by the test suite) into a Markdown evidence
appendix. Usage: python3 validation/report.py > build/validation-evidence.md"""
import json
import sys
from pathlib import Path

RES = Path(__file__).resolve().parents[1] / "build" / "validation-results"


def load(name):
    f = RES / f"{name}.jsonl"
    return [json.loads(x) for x in f.read_text().splitlines()] if f.exists() else []


def fmt(v):
    if isinstance(v, float):
        return f"{v:.3g}"
    if isinstance(v, list):
        return ", ".join(fmt(x) for x in v)
    return str(v)


def table(rows, cols):
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        out.append("| " + " | ".join(fmt(r.get(c, "")) for c in cols) + " |")
    return "\n".join(out)


def main():
    w = sys.stdout.write
    d = [r for r in load("test_differential") if "scenario" in r]
    w("## Differential: Fortran reference vs C++ port (SAC output)\n\n")
    w(f"{sum(r['files'] for r in d)} SAC files over {len(d)} scenarios; "
      f"{sum(r['bitwise_files'] for r in d)} bit-identical.\n\n")
    w(table(d, ["scenario", "files", "bitwise_files", "worst_rel_linf", "max_ulp", "headers_identical"]) + "\n\n")
    fm = [r for r in load("test_differential") if "fortran_fastmath_vs_strict" in r]
    if fm:
        w("Upstream default build (`-ffast-math`) vs strict reference, for scale: rel L-inf "
          f"{fmt(fm[0]['fortran_fastmath_vs_strict'])}; C++ vs fast-math build {fmt(fm[0]['cpp_vs_fortran_fastmath'])}.\n\n")
    for name, cols, title in [
        ("test_state", ["scenario", "checkpoints", "cell_values_compared", "bitwise"], "Full-wavefield checkpoints"),
        ("test_routines", ["cases", "mismatches"], "Routine-level probes"),
        ("test_determinism", ["test", "threads", "files", "bitwise"], "Thread-count determinism"),
        ("test_edge_cases", ["test", "fortran_rc", "cpp_rc", "rc", "stderr"], "Edge cases"),
        ("test_physics", ["test", "misfit", "observed_order", "max_relative_boundary_artifact", "limit",
                          "q_input", "q_estimated"], "Physics (independent of Fortran)"),
    ]:
        rows = load(name)
        if rows:
            w(f"## {title}\n\n" + table(rows, cols) + "\n\n")


if __name__ == "__main__":
    main()
