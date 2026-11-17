"""Readers for the full-wavefield checkpoint dumps written by SWPC_DUMP_DIR (Fortran + C++)."""
from __future__ import annotations

from pathlib import Path

import numpy as np

MEDIUM_FIELDS = ("rho", "lam", "mu", "taup", "taus")
STATE_FIELDS = (("Vx", "<f4"), ("Vz", "<f4"), ("Sxx", "<f8"), ("Szz", "<f8"), ("Sxz", "<f8"))


def _read(path: Path, fields):
    raw = Path(path).read_bytes()
    kb, ke, ib, ie = np.frombuffer(raw, "<i4", 4, 0)
    nk, ni = ke - kb + 1, ie - ib + 1
    off = 16
    out = {"bounds": (int(kb), int(ke), int(ib), int(ie))}
    for name, dt in fields:
        n = nk * ni
        a = np.frombuffer(raw, dt, n, off).reshape(ni, nk).T  # Fortran order a(k, i)
        out[name] = a.copy()
        off += n * np.dtype(dt).itemsize
    if off != len(raw):
        raise ValueError(f"{path}: {len(raw) - off} trailing bytes")
    return out


def read_medium(path):
    return _read(path, [(f, "<f4") for f in MEDIUM_FIELDS])


def read_state(path):
    return _read(path, STATE_FIELDS)


def steps(d: Path) -> list[int]:
    return sorted(int(p.stem.split("_")[1]) for p in Path(d).glob("state_*.bin"))
