"""Scenario definitions: generates OpenSWPC-compatible input (.inf), source, station and layer files."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_LHM = """# depth  rho  vp  vs  qp  qs
  0.0   2.3  4.0  2.3  200  100
  3.0   2.5  5.0  2.9  400  200
  8.0   2.8  6.2  3.6  600  300
"""


@dataclass
class Scenario:
    name: str
    nx: int = 200
    nz: int = 120
    nt: int = 400
    dx: float = 0.1
    dz: float = 0.1
    dt: float = 0.008
    xbeg: float | None = None
    zbeg: float = -1.0
    na: int = 20
    abc_type: str = "pml"
    fullspace_mode: bool = False
    vmodel_type: str = "uni"
    vp0: float = 5.0
    vs0: float = 3.0
    rho0: float = 2.7
    qp0: float = 1.0e6
    qs0: float = 1.0e6
    topo0: float = 0.0
    is_ocean: bool = False
    munk_profile: bool = False
    stabilize_pml: bool = False
    earth_flattening: bool = False
    lhm: str = DEFAULT_LHM
    fq_min: float = 0.05
    fq_max: float = 5.0
    fq_ref: float = 1.0
    bf_mode: bool = False
    stf_format: str = "xym0ij"
    stftype: str = "kupper"
    # each source row is the tail of a line after x y z (tokens as in OpenSWPC files)
    sources: list[str] = field(default_factory=lambda: ["0.0 0.0 3.0  0.0 0.5  1.0e15  1.0 0.0 1.0 0.0 0.0 0.0"])
    stations: list[str] = field(default_factory=lambda: [
        "-5.0 0.0 0.0 ST01 fsb", "-2.0 0.0 1.5 ST02 dep", "2.5 0.0 0.0 ST03 fsb", "6.0 0.0 5.0 ST04 dep"])
    sw_wav_v: bool = True
    sw_wav_u: bool = True
    sw_wav_stress: bool = True
    sw_wav_strain: bool = True
    ntdec_w: int = 2
    extra: dict[str, str] = field(default_factory=dict)

    def write(self, d: Path) -> Path:
        d.mkdir(parents=True, exist_ok=True)
        (d / "source.dat").write_text("# generated\n" + "\n".join(self.sources) + "\n")
        (d / "stloc.xy").write_text("# generated\n" + "\n".join(self.stations) + "\n")
        (d / "lhm.dat").write_text(self.lhm)
        xbeg = self.xbeg if self.xbeg is not None else -self.nx * self.dx / 2
        fl = lambda b: ".true." if b else ".false."  # noqa: E731
        kv = {
            "title": f"'{self.name}'", "odir": "'./out'", "ntdec_r": "100000", "nproc_x": "1",
            "nx": self.nx, "nz": self.nz, "nt": self.nt, "dx": self.dx, "dz": self.dz, "dt": self.dt,
            "xbeg": xbeg, "zbeg": self.zbeg, "tbeg": 0.0, "clon": 139.7604, "clat": 35.7182, "phi": 0.0,
            "fq_min": self.fq_min, "fq_max": self.fq_max, "fq_ref": self.fq_ref,
            "snp_format": "'native'", "xz_ps%sw": ".false.", "fs_ps%sw": ".false.", "ob_ps%sw": ".false.",
            "xz_v%sw": ".false.", "fs_v%sw": ".false.", "ob_v%sw": ".false.", "xz_u%sw": ".false.",
            "fs_u%sw": ".false.", "ob_u%sw": ".false.",
            "sw_wav_v": fl(self.sw_wav_v), "sw_wav_u": fl(self.sw_wav_u), "sw_wav_stress": fl(self.sw_wav_stress),
            "sw_wav_strain": fl(self.sw_wav_strain), "ntdec_w": self.ntdec_w, "st_format": "'xy'",
            "fn_stloc": "'./stloc.xy'", "wav_format": "'sac'", "ntdec_w_prg": 0,
            "stf_format": f"'{self.stf_format}'", "stftype": f"'{self.stftype}'", "fn_stf": "'./source.dat'",
            "sdep_fit": "'asis'", "bf_mode": fl(self.bf_mode), "pw_mode": ".false.",
            "abc_type": f"'{self.abc_type}'", "na": self.na, "stabilize_pml": fl(self.stabilize_pml),
            "fullspace_mode": fl(self.fullspace_mode), "vmodel_type": f"'{self.vmodel_type}'",
            "is_ocean": fl(self.is_ocean), "topo_flatten": ".false.", "munk_profile": fl(self.munk_profile),
            "earth_flattening": fl(self.earth_flattening),
            "vp0": self.vp0, "vs0": self.vs0, "rho0": self.rho0, "qp0": self.qp0, "qs0": self.qs0,
            "topo0": self.topo0, "fn_lhm": "'./lhm.dat'", "green_mode": ".false.", "stopwatch_mode": ".false.",
            "benchmark_mode": ".false.", "ipad": 0, "kpad": 0, "strict_mode": ".false.",
        }
        kv.update(self.extra)
        lines = [f"  {k:<18} = {v}" for k, v in kv.items()]
        inf = d / f"{self.name}.inf"
        inf.write_text("\n".join(lines) + "\n")
        return inf
