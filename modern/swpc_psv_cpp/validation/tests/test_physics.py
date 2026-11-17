"""Physics checks of the port that do not rely on the Fortran reference."""
import numpy as np
import pytest
from scipy.special import hankel2

from swpcval import runner
from swpcval.scenario import Scenario

VP, VS, RHO = 5.0, 3.0, 2.7
TS, TR = 0.1, 0.3


def kupper_moment(t):
    """Time integral of the Kupper moment-rate function (unit final moment)."""
    x = np.clip((t - TS) / TR, 0.0, 1.0)
    th = np.pi * x
    # integral of (3 pi / 4 tr) sin^3 -> (3/4)(2/3 - cos + cos^3/3)
    return 0.75 * (2.0 / 3.0 - np.cos(th) + np.cos(th) ** 3 / 3.0)


def analytic_radial_velocity(r, t, alpha=VP, nfft=1 << 15, dtf=0.001):
    """Radial particle velocity of a 2-D isotropic (explosion) line source in a full space,
    up to a constant factor: v_r(w) = -i k Mdot(w) H1^(2)(k r) / (4 rho alpha^2), k = w / alpha
    (time dependence exp(+i w t), so H^(2) is the outgoing wave)."""
    tt = np.arange(nfft) * dtf
    w = 2 * np.pi * np.fft.rfftfreq(nfft, dtf)
    mr = np.fft.rfft(np.gradient(kupper_moment(tt), dtf))  # moment rate, compact support
    k = w / alpha
    with np.errstate(invalid="ignore", divide="ignore"):
        g = np.where(w > 0, -1j * k * hankel2(1, k * r) / (4 * RHO * alpha ** 2), 0)  # v_r per unit moment rate
    v = np.fft.irfft(mr * g, nfft)
    return np.interp(t, tt, v)


def explosion(name, dx, nt_phys=2.2, stations=None, **kw):
    n = int(round(14.0 / dx))
    dt = 0.08 * dx
    st = stations or [f"{x:.3f} 0.0 {z:.3f} R{i:02d} dep" for i, (x, z) in enumerate(
        [(2.0, 0.0), (0.0, 3.0), (-2.5, -2.5), (3.5, 1.0), (-1.0, -4.0)])]
    return Scenario(name, nx=n, nz=n, dx=dx, dz=dx, dt=dt, nt=int(round(nt_phys / dt)), xbeg=-7.0, zbeg=-7.0,
                    fullspace_mode=True, vp0=VP, vs0=VS, rho0=RHO, sw_wav_u=False, sw_wav_stress=False,
                    sw_wav_strain=False, ntdec_w=1,
                    sources=[f"0.0 0.0 0.0  {TS} {TR}  1.0e15  1.0 0.0 1.0 0.0 0.0 0.0"], stations=st, **kw)


def misfit_vs_analytic(out, sc):
    """Single global amplitude scale fitted over all stations/components; returns (misfit, scale)."""
    num, ana = [], []
    pos = {}
    for line in sc.stations:
        x, _, z, nm, _ = line.split()
        pos[nm] = (float(x), float(z))
    for key, s in out.items():
        nm, cmp = key.split(".")
        x, z = pos[nm]
        r = np.hypot(x, z)
        # velocities live on the half time level of the leap-frog scheme
        vr = analytic_radial_velocity(r, s.time + sc.dt / 2)
        a = vr * x / r if cmp == "Vx" else -vr * z / r  # SAC Vz is positive upward, model z is down
        num.append(s.data.astype(float))
        ana.append(a)
    n, a = np.concatenate(num), np.concatenate(ana)
    scale = float(n @ a / (a @ a))
    return float(np.linalg.norm(n - scale * a) / np.linalg.norm(n)), scale


@pytest.mark.physics
def test_explosion_matches_analytic_solution(workdir, cpp_bin, record):
    sc = explosion("analytic", 0.025)
    runner.run("cpp", sc, workdir, threads=8)
    m, scale = misfit_vs_analytic(runner.outputs(workdir), sc)
    record(dx=0.025, misfit=m, fitted_scale=scale, stations=len(sc.stations))
    assert m < 0.03, f"waveform misfit vs analytic solution {m:.3%}"


@pytest.mark.physics
@pytest.mark.slow
def test_grid_convergence(workdir, cpp_bin, record):
    errs = []
    dxs = [0.1, 0.05, 0.025]
    for dx in dxs:
        sc = explosion(f"conv{int(dx * 1000)}", dx)
        runner.run("cpp", sc, workdir / f"dx{dx}", threads=8)
        errs.append(misfit_vs_analytic(runner.outputs(workdir / f"dx{dx}"), sc)[0])
    orders = [float(np.log2(errs[i] / errs[i + 1])) for i in range(len(errs) - 1)]
    record(dx=dxs, misfit=errs, observed_order=orders)
    # Error must fall monotonically with grid refinement. Below ~0.6 % the misfit stops being set by
    # the grid and is set by the point-source / analytic comparison floor, so only the coarse pair
    # is held to an order requirement.
    assert errs[0] > errs[1] > errs[2]
    assert orders[0] > 1.2, f"observed convergence orders {orders}"


@pytest.mark.physics
@pytest.mark.parametrize("abc,limit", [("pml", 0.02), ("cerjan", 0.15)])
def test_absorbing_boundary_reflection(abc, limit, workdir, cpp_bin, record):
    """Compare a small domain against a domain large enough that no boundary signal returns."""
    st = ["3.0 0.0 0.0 RA dep", "2.0 0.0 2.0 RB dep", "-3.0 0.0 -1.5 RC dep"]
    kw = dict(dx=0.05, dz=0.05, dt=0.004, nt=1000, fullspace_mode=True, sw_wav_u=False, sw_wav_stress=False,
              sw_wav_strain=False, ntdec_w=1, abc_type=abc, stations=st,
              sources=[f"0.0 0.0 0.0  {TS} {TR}  1.0e15  1.0 0.0 0.5 0.0 0.4 0.0"])
    small = Scenario("small", nx=200, nz=200, xbeg=-5.0, zbeg=-5.0, **kw)
    big = Scenario("big", nx=640, nz=640, xbeg=-16.0, zbeg=-16.0, **kw)
    runner.run("cpp", small, workdir / "s", threads=8)
    runner.run("cpp", big, workdir / "b", threads=8)
    so, bo = runner.outputs(workdir / "s"), runner.outputs(workdir / "b")
    worst = max(float(np.max(np.abs(so[k].data - bo[k].data)) / np.max(np.abs(bo[k].data))) for k in so)
    record(abc=abc, na=20, max_relative_boundary_artifact=worst, limit=limit)
    assert worst < limit


@pytest.mark.physics
def test_reciprocity_in_layered_medium(workdir, cpp_bin, record):
    """G_xz(B; A) = G_zx(A; B) for point forces in a heterogeneous, attenuating layered model."""
    A, B = (-1.5, 1.0), (2.0, 3.5)
    common = dict(vmodel_type="lhm", bf_mode=True, stf_format="xy", nx=300, nz=200, dx=0.05, dz=0.05, dt=0.004,
                  nt=900, zbeg=-1.0, sw_wav_u=False, sw_wav_stress=False, sw_wav_strain=False, ntdec_w=1)
    # body-force rows: x y z ts tr fx fy fz
    s1 = Scenario("recip1", sources=[f"{A[0]} 0.0 {A[1]}  0.0 0.3  0.0 0.0 1.0e15"],
                  stations=[f"{B[0]} 0.0 {B[1]} RB dep"], **common)
    s2 = Scenario("recip2", sources=[f"{B[0]} 0.0 {B[1]}  0.0 0.3  1.0e15 0.0 0.0"],
                  stations=[f"{A[0]} 0.0 {A[1]} RA dep"], **common)
    runner.run("cpp", s1, workdir / "a", threads=8)
    runner.run("cpp", s2, workdir / "b", threads=8)
    gxz = runner.outputs(workdir / "a")["RB.Vx"].data.astype(float)
    gzx = -runner.outputs(workdir / "b")["RA.Vz"].data.astype(float)
    err = float(np.linalg.norm(gxz - gzx) / np.linalg.norm(gxz))
    record(misfit=err)
    assert err < 0.01


@pytest.mark.physics
def test_constant_q_attenuation(workdir, cpp_bin, record):
    """Spectral ratio between a Q=50 and an elastic run at two distances recovers Q."""
    Q = 50.0
    st = ["2.0 0.0 0.0 N1 dep", "5.0 0.0 0.0 N2 dep"]
    kw = dict(nx=440, nz=280, dx=0.04, dz=0.04, dt=0.004, nt=900, xbeg=-3.0, zbeg=-5.6, fullspace_mode=True,
              stations=st, sw_wav_u=False, sw_wav_stress=False, sw_wav_strain=False, ntdec_w=1,
              fq_min=0.1, fq_max=10.0, fq_ref=1.0, sources=["0.0 0.0 0.0  0.0 0.2  1.0e15  1.0 0.0 1.0 0.0 0.0 0.0"])
    runner.run("cpp", Scenario("q", qp0=Q, qs0=Q, **kw), workdir / "q", threads=8)
    runner.run("cpp", Scenario("el", qp0=1e6, qs0=1e6, **kw), workdir / "e", threads=8)
    q, e = runner.outputs(workdir / "q"), runner.outputs(workdir / "e")
    dt = q["N1.Vx"].f("delta")
    f = np.fft.rfftfreq(4096, dt)
    spec = lambda s: np.abs(np.fft.rfft(s.data.astype(float), 4096))  # noqa: E731
    lr = np.log((spec(q["N2.Vx"]) / spec(e["N2.Vx"])) / (spec(q["N1.Vx"]) / spec(e["N1.Vx"])))
    band = (f > 0.5) & (f < 4.0)
    slope = np.polyfit(f[band], lr[band], 1)[0]
    q_est = -np.pi * (5.0 - 2.0) / VP / slope
    record(q_input=Q, q_estimated=float(q_est), tolerance=0.25)
    # The solver approximates constant Q with NM = 3 relaxation mechanisms; the estimate is
    # band-dependent (53-62 across 0.5-6 Hz bands for Q = 50), hence the 25 % tolerance.
    assert abs(q_est - Q) / Q < 0.25
