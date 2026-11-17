"""Input handling, stability guard and scope boundaries: both codes must agree on what is
accepted and what is rejected; the port must fail loudly (non-zero exit) outside its scope."""
import pytest

from swpcval import runner
from swpcval.scenario import Scenario


def _rc(engine, sc, d):
    return runner.run(engine, sc, d, check=False).returncode


@pytest.mark.differential
def test_cfl_violation_rejected_by_both(workdir, cpp_bin, fortran_ref, record):
    sc = Scenario("cfl_bad", dt=0.02, nt=50)  # c = vp*dt*(7/6)*sqrt(2)/dx > 1
    rc_f = _rc("fortran", sc, workdir / "f")
    rc_c = _rc("cpp", sc, workdir / "c")
    p = runner.run("cpp", sc, workdir / "c", check=False)
    record(fortran_rc=rc_f, cpp_rc=rc_c)
    assert rc_f != 0 and rc_c != 0
    assert "CFL" in p.stderr
    assert not list((workdir / "c" / "out" / "wav").glob("*.sac"))


@pytest.mark.differential
def test_cfl_marginal_accepted_by_both(workdir, cpp_bin, fortran_ref):
    sc = Scenario("cfl_ok", dt=0.0105, nt=50)  # c just below 1
    assert _rc("fortran", sc, workdir / "f") == 0
    assert _rc("cpp", sc, workdir / "c") == 0


@pytest.mark.parametrize("extra,msg", [
    ({"vmodel_type": "'grd'"}, "vmodel_type"),
    ({"vmodel_type": "'lhm_rmed'"}, "vmodel_type"),
    ({"pw_mode": ".true."}, "pw_mode"),
    ({"wav_format": "'csf'"}, "wav_format"),
    ({"stf_format": "'bogus'"}, "stf_format"),
    ({"fn_stf": "'./missing.dat'"}, "cannot open"),
])
def test_out_of_scope_inputs_fail_loudly(extra, msg, workdir, cpp_bin, record):
    p = runner.run("cpp", Scenario("scope", nt=10, extra=extra), workdir, check=False)
    record(extra=extra, rc=p.returncode, stderr=p.stderr.strip()[:200])
    assert p.returncode != 0 and msg in p.stderr


def test_missing_parameter_file(workdir, cpp_bin):
    import subprocess
    p = subprocess.run([str(cpp_bin), "-i", str(workdir / "nope.inf")], capture_output=True, text=True)
    assert p.returncode != 0


def test_no_station_file_means_no_output(workdir, cpp_bin):
    sc = Scenario("nost", nt=20, extra={"fn_stloc": "'./absent.xy'"})
    runner.run("cpp", sc, workdir)
    assert not list((workdir / "out" / "wav").glob("*.sac"))


@pytest.mark.differential
def test_parser_semantics_match_fortran(workdir, cpp_bin, fortran_ref):
    """Duplicate keys (first wins), comments, Fortran D-exponents, logical spellings, defaults."""
    sc = Scenario("parse", nt=120)
    inf = sc.write(workdir / "f")
    text = inf.read_text().replace("  nt                 = 120", "  nt                 = 120   ! trailing comment\n  nt = 999")
    text = text.replace("  sw_wav_v           = .true.", "  sw_wav_v           = T")
    text = text.replace("  dt                 = 0.008", "  dt                 = 8.0D-3")
    text = "# leading comment line\n! another\n" + text.replace("  ntdec_w            = 2\n", "")  # default 10
    import shutil
    from swpcval.sac import read_sac
    outs = {}
    for eng in ("fortran", "cpp"):
        d = workdir / eng
        if d.exists():
            shutil.rmtree(d)
        sc.write(d).write_text(text)
        import os
        import subprocess
        env = dict(os.environ, OMP_NUM_THREADS="1", HWLOC_COMPONENTS="no_os,stop", MALLOC_MMAP_THRESHOLD_="0")
        cmd = ([str(cpp_bin), "-i", "parse.inf", "-q"] if eng == "cpp" else
               ["mpirun", "--oversubscribe", "--allow-run-as-root", "-np", "1", str(fortran_ref), "-i", "parse.inf"])
        subprocess.run(cmd, cwd=d, env=env, check=True, capture_output=True)
        outs[eng] = {f.name: read_sac(f) for f in (d / "out" / "wav").glob("*.sac")}
    assert sorted(outs["fortran"]) == sorted(outs["cpp"])
    for k, f in outs["fortran"].items():
        c = outs["cpp"][k]
        assert f.data.tobytes() == c.data.tobytes()
        assert f.i("npts") == c.i("npts") == 12  # nt=120 (first wins), ntdec_w default 10
