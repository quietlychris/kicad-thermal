#!/usr/bin/env python3
"""Run kicad_fdtd (Python) and kicad-thermal (Rust) on the same geometry and compare.

    python scripts/compare_live.py board_geometry.json [--kicad-fdtd PATH] \\
        -- --power U3=1.8 R12=0.4 --h 5 --tamb 45 --cell 0.25 \\
           --net /12V --from U3.5 --to J7.1 --current 5

Everything after `--` goes to both tools unchanged (the `kicad-fdtd thermal` flags).
Use this on real exported boards, which the synthetic test fixtures can't fully stand
in for. Prints wall times, the worst field differences normalized like the regression
suite (temperature by the max rise, voltage by the drop, current density by Jmax), and
whether the report text matches. Exits non-zero if any difference exceeds --tol.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
import time

import numpy as np

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def load(path):
    with open(path) as f:
        return json.load(f)


def worst(a, b, layers, mask_value=None):
    w = 0.0
    for l in layers:
        x, y = np.array(a[l], float), np.array(b[l], float)
        if mask_value is not None:
            if not np.array_equal(x == mask_value, y == mask_value):
                return float("inf")
            keep = x != mask_value
            x, y = x[keep], y[keep]
        if x.size:
            w = max(w, float(np.abs(x - y).max()))
    return w


def main():
    if "--" in sys.argv:
        i = sys.argv.index("--")
        own, passthru = sys.argv[1:i], sys.argv[i + 1:]
    else:
        own, passthru = sys.argv[1:], []
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("geometry")
    ap.add_argument("--kicad-fdtd", default=os.environ.get("KICAD_FDTD") or os.path.join(ROOT, "third_party", "kicad_fdtd"),
                    help="kicad_fdtd checkout with the AMG patch (scripts/fetch_kicad_fdtd.sh; default third_party/kicad_fdtd)")
    ap.add_argument("--rust-bin", default=os.path.join(ROOT, "target", "release", "kicad-thermal"))
    ap.add_argument("--python-precond", default="jacobi", choices=["jacobi", "amg", "auto"])
    ap.add_argument("--rust-precond", default="amg", choices=["jacobi", "amg", "auto"])
    ap.add_argument("--tol", type=float, default=1e-6)
    a = ap.parse_args(own)
    if not os.path.exists(a.rust_bin):
        sys.exit("no Rust binary at %s: run `cargo build --release` first" % a.rust_bin)

    tmp = tempfile.mkdtemp(prefix="kt_compare_")
    py_out, rs_out = os.path.join(tmp, "py"), os.path.join(tmp, "rs")
    env = dict(os.environ)
    if a.kicad_fdtd:
        env["PYTHONPATH"] = a.kicad_fdtd + os.pathsep + env.get("PYTHONPATH", "")
    py_cmd = [sys.executable, "-m", "kicad_fdtd.thermal", a.geometry, *passthru, "--cpu", "--precond", a.python_precond, "--out", py_out]
    rs_cmd = [a.rust_bin, a.geometry, *passthru, "--precond", a.rust_precond, "--out", rs_out, "-v"]
    times = {}
    for name, cmd in (("python-" + a.python_precond, py_cmd), ("rust-" + a.rust_precond, rs_cmd)):
        t = time.time()
        p = subprocess.run(cmd, env=env, capture_output=True, text=True)
        times[name] = time.time() - t
        if p.returncode != 0:
            sys.exit("%s failed:\n%s%s" % (name, p.stdout, p.stderr))
        if name.startswith("rust") and p.stderr.strip():
            print(p.stderr.strip())
    py = load(os.path.join(py_out, os.listdir(py_out)[0]))
    rs = load(os.path.join(rs_out, os.listdir(rs_out)[0]))

    layers = py["layers"]
    rise = max(abs(py["Tmax"] - py["tamb"]), 1e-12)
    rows = [("temperature (of max rise)", worst(py["T"], rs["T"], layers) / rise)]
    if "dc" in py:
        drop = max(abs(py["dc"]["drop_V"]), 1e-30)
        rows.append(("DC voltage (of drop)", worst(py["dc"]["V"], rs["dc"]["V"], layers, mask_value=-1.0) / drop))
        rows.append(("current density (of Jmax)", worst(py["dc"]["J"], rs["dc"]["J"], layers) / max(py["dc"]["Jmax"], 1e-30)))
    for name, t in times.items():
        print("%-16s %7.1f s wall" % (name, t))
    ok = True
    for name, w in rows:
        flag = "ok" if w <= a.tol else "EXCEEDS --tol %g" % a.tol
        ok &= w <= a.tol
        print("max difference, %-27s %.2e  %s" % (name, w, flag))
    same = py["report"] == rs["report"]
    print("report text identical: %s" % same)
    if not same:
        for lp, lr in zip(py["report"].splitlines(), rs["report"].splitlines()):
            if lp != lr:
                print("  python: %s\n  rust:   %s" % (lp, lr))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
