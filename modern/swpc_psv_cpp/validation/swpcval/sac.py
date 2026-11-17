"""Minimal binary SAC (v6, native little-endian) reader used by the validation suite."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

F_NAMES = {1: "delta", 2: "depmin", 3: "depmax", 6: "b", 7: "e", 8: "o", 32: "stla", 33: "stlo", 35: "stdp",
           36: "evla", 37: "evlo", 39: "evdp", 40: "mag", 41: "user0", 43: "user2", 45: "user4", 47: "user6",
           48: "user7", 49: "user8", 51: "dist", 52: "az", 53: "baz", 57: "depmen", 58: "cmpaz", 59: "cmpinc"}
I_NAMES = {71: "nzyear", 72: "nzjday", 73: "nzhour", 74: "nzmin", 75: "nzsec", 76: "nzmsec", 77: "nvhdr",
           80: "npts", 86: "iftype", 87: "idep"}
# Header words that legitimately differ between two runs (wall-clock execution date).
TIME_WORDS = {71, 72, 73, 74, 75}


@dataclass
class Sac:
    floats: np.ndarray    # 70 float32
    ints: np.ndarray      # 40 int32 (words 71..110 incl. logicals)
    chars: bytes          # 192 bytes
    data: np.ndarray      # float32

    def f(self, name: str) -> float:
        return float(self.floats[{v: k for k, v in F_NAMES.items()}[name] - 1])

    def i(self, name: str) -> int:
        return int(self.ints[{v: k for k, v in I_NAMES.items()}[name] - 71])

    @property
    def kstnm(self) -> str:
        return self.chars[0:8].decode().strip()

    @property
    def kcmpnm(self) -> str:
        return self.chars[160:168].decode().strip()

    @property
    def time(self) -> np.ndarray:
        return self.f("b") + self.f("delta") * np.arange(self.data.size)


def read_sac(path: str | Path) -> Sac:
    raw = Path(path).read_bytes()
    if len(raw) < 632:
        raise ValueError(f"{path}: truncated SAC file ({len(raw)} bytes)")
    floats = np.frombuffer(raw, "<f4", 70, 0).copy()
    ints = np.frombuffer(raw, "<i4", 40, 280).copy()
    chars = raw[440:632]
    npts = int(ints[80 - 71])
    data = np.frombuffer(raw, "<f4", npts, 632).copy()
    if 632 + 4 * npts != len(raw):
        raise ValueError(f"{path}: size {len(raw)} inconsistent with npts={npts}")
    return Sac(floats, ints, chars, data)


def header_diff(a: Sac, b: Sac, ignore_time: bool = True) -> list[str]:
    """Return a list of human-readable header mismatches (bitwise for float/int words)."""
    out = []
    for w in range(70):
        if a.floats[w].tobytes() != b.floats[w].tobytes():
            out.append(f"word {w + 1} ({F_NAMES.get(w + 1, '?')}): {a.floats[w]!r} != {b.floats[w]!r}")
    for w in range(40):
        if ignore_time and (w + 71) in TIME_WORDS:
            continue
        if a.ints[w] != b.ints[w]:
            out.append(f"word {w + 71} ({I_NAMES.get(w + 71, '?')}): {a.ints[w]} != {b.ints[w]}")
    if a.chars != b.chars:
        out.append(f"character block differs: {a.chars!r} != {b.chars!r}")
    return out
