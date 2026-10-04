#!/usr/bin/env python3
"""Generate the Python reference results the Rust regression suite compares against.

    scripts/fetch_kicad_fdtd.sh && python scripts/gen_reference.py [--kicad-fdtd PATH]

Runs kicad_fdtd's thermal module (commit 74e21eb plus the AMG preconditioner patch, so
both `precond="jacobi"` -- the original solver -- and `precond="amg"` are available) on
every case in tests/fixtures/cases.json and writes tests/reference/<case>/: the grid
materials (kx, kz, via copper area, copper fractions, inside-outline mask) and, for each
preconditioner, the full temperature, DC voltage and current-density fields plus the
scalar results, report text and CG iteration counts. Arrays are little-endian float64,
gzip-compressed with a fixed timestamp so regenerated files diff cleanly.

Needs numpy, scipy, shapely and pyamg. Only rerun this when the Python reference itself
changes; the Rust suite treats these files as ground truth.
"""
import argparse
import gzip
import json
import os
import platform
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
FIX = os.path.join(ROOT, "tests", "fixtures")
REF = os.path.join(ROOT, "tests", "reference")


def save(dirpath, name, arr):
    path = os.path.join(dirpath, name + ".f64.gz")
    data = np.ascontiguousarray(arr, dtype="<f8").tobytes()
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as f:
        f.write(data)
    return os.path.basename(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kicad-fdtd", default=os.environ.get("KICAD_FDTD") or os.path.join(ROOT, "third_party", "kicad_fdtd"),
                    help="kicad_fdtd checkout with the AMG patch (scripts/fetch_kicad_fdtd.sh; default third_party/kicad_fdtd)")
    a = ap.parse_args()
    if a.kicad_fdtd:
        sys.path.insert(0, a.kicad_fdtd)
    from kicad_fdtd import geometry, thermal
    import scipy, shapely, pyamg
    import scipy.sparse.linalg as sla
    if "precond" not in thermal._cg.__code__.co_varnames:
        sys.exit("this kicad_fdtd has no precond option: apply kicad_fdtd-amg-preconditioner.patch first")

    # count CG iterations through scipy's callback
    count = [0]
    orig_cg = sla.cg
    def counting_cg(*args, **kw):
        kw["callback"] = lambda xk: count.__setitem__(0, count[0] + 1)
        return orig_cg(*args, **kw)
    sla.cg = counting_cg

    commit = None
    if a.kicad_fdtd:
        try:
            commit = subprocess.run(["git", "-C", a.kicad_fdtd, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip() or None
        except OSError:
            pass
    generator = {"kicad_fdtd_commit": commit, "python": platform.python_version(), "numpy": np.__version__,
                 "scipy": scipy.__version__, "shapely": shapely.__version__, "pyamg": pyamg.__version__}

    cases = json.load(open(os.path.join(FIX, "cases.json")))
    for case in cases:
        name = case["name"]
        out = os.path.join(REF, name)
        os.makedirs(out, exist_ok=True)
        bd = geometry.Board(os.path.join(FIX, name + ".json"))
        g = thermal.Grid(bd, cell=case["cell"])
        meta = {"case": case, "generator": generator, "grid": {"nx": g.nx, "ny": g.ny, "nz": g.nz}, "arrays": {}, "results": {}}
        meta["arrays"]["kx"] = save(out, "kx", g.kx)
        meta["arrays"]["kz"] = save(out, "kz", g.kz)
        meta["arrays"]["via_area"] = save(out, "via_area", g.via_copper_area())
        meta["arrays"]["inside"] = save(out, "inside", g.inside_outline().astype(float))
        for l, layer in enumerate(bd.layers):
            meta["arrays"]["cu_frac_%d" % l] = save(out, "cu_frac_%d" % l, g.copper_fraction(layer))
        powers = {r: w for r, w in case["powers"]}
        for pc in ("jacobi", "amg"):
            res = {"arrays": {}}
            dc = None
            if case["dc"]:
                d = case["dc"]
                count[0] = 0
                dc = thermal.solve_dc(g, d["net"], tuple(d["from"]), tuple(d["to"]), d["current"], use_gpu=False, precond=pc)
                res["dc"] = {"drop_V": dc["drop_V"], "R_ohm": dc["R_ohm"], "joule_W": dc["joule_W"], "Jmax": dc["Jmax"],
                             "cells": dc["cells"], "iterations": count[0]}
                res["arrays"]["V"] = save(out, pc + "_V", dc["V"])
                res["arrays"]["J"] = save(out, pc + "_J", dc["J"])
            count[0] = 0
            th = thermal.solve_thermal(g, powers, h=case["h"], tamb=case["tamb"], extra_q=dc["joule"] if dc else None,
                                       use_gpu=False, precond=pc)
            res["thermal"] = {"Tmax": th["Tmax"], "P_total": th["P_total"], "iterations": count[0]}
            res["arrays"]["T"] = save(out, pc + "_T", th["T"])
            res["report"] = thermal.report(g, th, dc)
            meta["results"][pc] = res
        with open(os.path.join(out, "meta.json"), "w") as f:
            json.dump(meta, f, indent=1)
        print("%-15s %4d x %4d x %2d   jacobi %s it, amg %s it" % (name, g.nx, g.ny, g.nz,
              meta["results"]["jacobi"]["thermal"]["iterations"], meta["results"]["amg"]["thermal"]["iterations"]))


if __name__ == "__main__":
    main()
