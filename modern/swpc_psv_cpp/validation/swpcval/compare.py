"""Waveform comparison metrics used by the differential tests."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Metrics:
    max_abs: float        # max |a - b|
    peak: float           # max |a|
    rel_linf: float       # max|a-b| / max|a|
    rel_l2: float         # ||a-b||_2 / ||a||_2
    max_ulp: int          # max float32 ULP distance
    bitwise: bool

    def row(self) -> dict:
        return self.__dict__.copy()


def ulp_distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    ia = a.astype("<f4").view("<i4").astype(np.int64)
    ib = b.astype("<f4").view("<i4").astype(np.int64)
    ia = np.where(ia < 0, np.int64(-2**31) - ia, ia)
    ib = np.where(ib < 0, np.int64(-2**31) - ib, ib)
    return np.abs(ia - ib)


def compare(a: np.ndarray, b: np.ndarray) -> Metrics:
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch {a.shape} vs {b.shape}")
    d = np.abs(a.astype(np.float64) - b.astype(np.float64))
    peak = float(np.max(np.abs(a))) if a.size else 0.0
    na = float(np.linalg.norm(a.astype(np.float64)))
    return Metrics(
        max_abs=float(d.max()) if d.size else 0.0,
        peak=peak,
        rel_linf=float(d.max() / peak) if peak > 0 else float(d.max() if d.size else 0.0),
        rel_l2=float(np.linalg.norm(d) / na) if na > 0 else float(np.linalg.norm(d)),
        max_ulp=int(ulp_distance(a, b).max()) if a.size else 0,
        bitwise=bool(a.tobytes() == b.tobytes()),
    )
