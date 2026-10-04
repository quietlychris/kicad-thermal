#!/usr/bin/env python3
"""cmp_grid.py NAME CELL: compare the grid materials (kx, kz, via copper area, copper fractions, outline mask) of a
benchmark board between kicad-thermal (examples/dump_grid) and the Python original, cell by cell. A debugging aid for
model changes; the regression suite does the same check on the test fixtures.
Needs `cargo build --release --examples`, bench/boards.py, and scripts/fetch_kicad_fdtd.sh (or $KICAD_FDTD)."""
import os, subprocess, sys, tempfile
import numpy as np
KT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.environ.get("KICAD_FDTD") or os.path.join(KT, "third_party", "kicad_fdtd"))
from kicad_fdtd import geometry, thermal  # noqa: E402
BOARDS = os.environ.get("BENCH_BOARDS") or os.path.join(KT, "bench", "boards")
name, cell = sys.argv[1], float(sys.argv[2])
geo, dump = os.path.join(BOARDS, name + ".json"), tempfile.mkdtemp(prefix="kt_dump_")
subprocess.run([os.path.join(KT, "target", "release", "examples", "dump_grid"), geo, str(cell), dump], check=True, capture_output=True)
bd = geometry.Board(geo); g = thermal.Grid(bd, cell=cell)
sh = (g.nx, g.ny, g.nz); r = lambda n: np.fromfile(os.path.join(dump, n + ".f64"))
out = []
for nm, py in (("kx", g.kx), ("kz", g.kz), ("via", g.via_copper_area())):
    out.append("%s diff cells %d" % (nm, (np.abs(r(nm).reshape(sh) - py) > 0).sum()))
for l, layer in enumerate(bd.layers):
    out.append("%s cufrac diff %d" % (layer, (r("cufrac_%d" % l).reshape(g.nx, g.ny) != g.copper_fraction(layer)).sum()))
out.append("inside diff %d" % (r("inside").reshape(g.nx, g.ny).astype(bool) != g.inside_outline()).sum())
print(name, cell, " | ".join(out))
