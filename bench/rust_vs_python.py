#!/usr/bin/env python3
"""Rust (release) vs Python benchmark + agreement check on the README benchmark boards.

    scripts/fetch_kicad_fdtd.sh && cargo build --release --locked && python3 bench/boards.py && python3 bench/rust_vs_python.py

Each tool runs once per config, sequentially; times are whole-process wall clock. Every run is
diffed against Python Jacobi (the baseline) using compare_live.py's metrics and the 1e-6
tolerance; exits non-zero if any Rust run exceeds it or its report text differs.
Python runs with single-threaded BLAS; Rust at each RAYON_NUM_THREADS in $THREADS.
Boards come from $BENCH_BOARDS (default bench/boards, where boards.py writes them); the Python original from
$KICAD_FDTD (default third_party/kicad_fdtd, see scripts/fetch_kicad_fdtd.sh);
outputs and results.json go to $BENCH_OUT (default a fresh temp directory).
"""
import json, os, subprocess, sys, tempfile, time
KT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
BOARD_DIR = os.environ.get("BENCH_BOARDS") or os.path.join(KT, "bench", "boards")
KFDTD = os.environ.get("KICAD_FDTD") or os.path.join(KT, "third_party", "kicad_fdtd")
S = os.environ.get("BENCH_OUT") or tempfile.mkdtemp(prefix="kt_bench_")
RUST = KT + "/target/release/kicad-thermal"
BOARDS = [("B_2L_compact", 0.25), ("D_2L_long_strip", 0.25), ("A_4L_power", 0.25),
          ("C_6L_heavy", 0.25), ("B_2L_compact", 0.1)]
THREADS = [int(t) for t in os.environ.get("THREADS", "1,4,16").split(",")]
sys.path.insert(0, KT + "/scripts")
from compare_live import load, worst  # noqa: E402
ONE = {k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}


def args_for(geo):
    g = json.load(open(geo))
    us = [f["ref"] for f in g["footprints"] if f["ref"].startswith("U")]
    pw = ["%s=%.2f" % (u, 0.15 + 0.05 * i) for i, u in enumerate(us)]
    return ["--power", *pw, "--h", "5", "--tamb", "45", "--net", "/PWR", "--from", "S.1", "--to", "K.1", "--current", "2"]


def run(cmd, env, out):
    t = time.time()
    p = subprocess.run(cmd, env=env, capture_output=True, text=True)
    dt = time.time() - t
    if p.returncode:
        sys.exit("FAILED %s\n%s%s" % (cmd, p.stdout, p.stderr))
    return dt, load(os.path.join(out, os.listdir(out)[0])), p.stderr


def diff(py, rs):
    L = py["layers"]
    d = {"T": worst(py["T"], rs["T"], L) / max(abs(py["Tmax"] - py["tamb"]), 1e-12)}
    if "dc" in py:
        d["V"] = worst(py["dc"]["V"], rs["dc"]["V"], L, mask_value=-1.0) / max(abs(py["dc"]["drop_V"]), 1e-30)
        d["J"] = worst(py["dc"]["J"], rs["dc"]["J"], L) / max(py["dc"]["Jmax"], 1e-30)
    d["report"] = py["report"] == rs["report"]
    return d


rows = []
for name, cell in BOARDS:
    geo = "%s/%s.json" % (BOARD_DIR, name)
    a = args_for(geo) + ["--cell", str(cell)]
    tag = "%s@%g" % (name, cell)
    res = {}
    for pc in ("jacobi", "amg"):
        out = "%s/out/%s/py-%s" % (S, tag, pc); os.makedirs(out, exist_ok=True)
        env = dict(os.environ, PYTHONPATH=KFDTD, **ONE)
        res["py-" + pc] = run([sys.executable, "-m", "kicad_fdtd.thermal", geo, *a, "--cpu", "--precond", pc, "--out", out], env, out)
        print("%-22s py-%-6s       %7.2f s" % (tag, pc, res["py-" + pc][0]), flush=True)
    for pc in ("jacobi", "amg"):
        for th in THREADS:
            k = "rs-%s-t%d" % (pc, th)
            out = "%s/out/%s/%s" % (S, tag, k); os.makedirs(out, exist_ok=True)
            env = dict(os.environ, RAYON_NUM_THREADS=str(th))
            res[k] = run([RUST, geo, *a, "--precond", pc, "--out", out, "-v"], env, out)
            print("%-22s %-14s %7.2f s" % (tag, k, res[k][0]), flush=True)
    base = res["py-jacobi"][1]
    for k, (dt, r, err) in res.items():
        d = diff(base, r) if k != "py-jacobi" else {}
        rows.append({"board": tag, "run": k, "wall_s": dt, **d, "stderr": err if k.startswith("rs") else ""})
json.dump(rows, open(S + "/results.json", "w"), indent=1)

print("\n%-22s %-14s %8s %9s %9s %9s %s" % ("board", "run", "wall s", "dT/rise", "dV/drop", "dJ/Jmax", "report"))
bad = 0
for r in rows:
    if "T" in r:
        ok = max(r["T"], r.get("V", 0), r.get("J", 0)) <= 1e-6 and r["report"]
        # py-amg is informational; the contract is Rust vs Python Jacobi
        bad += (not ok) and r["run"].startswith("rs")
        print("%-22s %-14s %8.2f %9.1e %9.1e %9.1e %s%s" % (r["board"], r["run"], r["wall_s"], r["T"], r.get("V", 0), r.get("J", 0),
              "same" if r["report"] else "DIFFERS", "" if ok else "  <-- FAIL"))
    else:
        print("%-22s %-14s %8.2f   (baseline)" % (r["board"], r["run"], r["wall_s"]))
print("\nRust runs outside 1e-6 or with differing report: %d\nresults: %s/results.json" % (bad, S))
sys.exit(1 if bad else 0)
