# kicad-thermal

Steady-state thermal and DC IR-drop analysis on a KiCad board's exact copper, in Rust.

kicad-thermal is a port of the thermal module of [kicad_fdtd](https://github.com/Jerome-Graves/kicad_fdtd)
by Jerome Graves (MIT), at commit `74e21eb`. The physical model is unchanged and the results match the
original (see [Validation](#validation)). The port adds an algebraic multigrid (AMG) preconditioner, which on
real boards makes the solve about 30 times faster.

## What it computes

- **Thermal:** finite-volume heat conduction on a whole-board grid. Every copper layer is one cell thick, with
  its conductivity scaled by its copper fraction (the exact polygons are sampled 3 × 3 per cell). FR4 is
  0.8 W/mK in plane and 0.3 W/mK through the board; via barrels add through-plane copper. Component power goes
  into the cells under each footprint's pads. Both faces lose heat by convection *h* to *T_amb*; the edges are
  adiabatic, with no radiation and no enclosure.
- **DC:** sheet conduction in one net's copper plus via barrels, with current in at one pad and out at another.
  It reports the IR drop, resistance, current density per cell and Joule heat. The Joule heat is added to the
  thermal sources.

The model's limits are kicad_fdtd's. Temperatures are a conservative ranking of hot spots, not junction
temperatures. Peak temperatures and the voltage drop in narrow necks depend on the cell size, so check that a
result stops changing as you refine the grid.

## Install

```bash
cargo install --path .          # or: cargo build --release  (binary in target/release/kicad-thermal)
```

## Use

kicad-thermal reads a geometry JSON exported from a `.kicad_pcb`. The exporter needs KiCad's Python module
`pcbnew`, so that step stays in Python. `scripts/export_geometry.py` writes it with KiCad 9's own API and needs
nothing else:

```bash
# 1. Export the board geometry (--fill refills zones first; --stackup takes the stackup from a JSON file)
python3 scripts/export_geometry.py board.kicad_pcb board_geometry.json --fill

# 2. Solve
kicad-thermal board_geometry.json \
    --power U3=1.8 R12=0.4 --h 5 --tamb 45 --cell 0.25 \
    --net /12V --from U3.5 --to J7.1 --current 5 -v
```

The flags match `kicad-fdtd thermal`, plus:

- `--precond auto|amg|jacobi` (default `auto`): `jacobi` is the original solver and `amg` is smoothed-aggregation
  multigrid. `auto` uses AMG for systems of at least 50,000 unknowns and Jacobi below that, where AMG setup
  costs more than it saves.
- `-v` prints the number of unknowns, iterations, setup and solve times, and the AMG hierarchy to stderr.
- `--cpu` is accepted for compatibility and ignored; kicad-thermal always runs on the CPU.

The report printed to stdout, and the `out/thermal_<board>.json` result, have the same format as kicad_fdtd's.
The library ([`Board::load`](src/geometry.rs) → [`Grid::new`](src/grid.rs) → [`solve_dc`](src/solve.rs) →
[`solve_thermal`](src/solve.rs)) mirrors kicad_fdtd's Python API: build one grid, solve many DC paths on it, and
pass their summed Joule heat to a single thermal solve (`extra_q`).

**The geometry format** is the one kicad_fdtd's exporter (`kicad-fdtd export`) writes, and `src/geometry.rs` reads
either. `scripts/export_geometry.py` adds two things: each via's pad as copper on its top and bottom layers (a
track ending inside a via pad but outside its drill is connected in KiCad, and looked open without it), and
`--fill`/`--stackup`. On a 4-layer, 200-part board its pads, tracks, zones, vias, pads list, outline and stackup were
identical to kicad_fdtd's export; only the added via pads differ.

**KiCad 9 and kicad_fdtd's exporter:** at `74e21eb` it calls `BOARD.GetBoardPolygonOutlines` with KiCad 8's
signature. `scripts/kicad9_shim/sitecustomize.py` adapts that one call at runtime when it is on `PYTHONPATH`, for
comparisons against the Python original.

## Validation

**Test suite** (`cargo test --release`, no Python needed):

- `tests/regression.rs` checks seven synthetic boards against reference results from the Python original
  (`tests/reference/`). They cover 2-, 4- and 6-layer stackups, heavy copper, overlapping plane cutouts, an
  L-shaped outline with air cells, a board with no stackup, a DC path through vias, through-hole and non-plated
  pads, and kicad_fdtd's two analytic cases. Each case checks four things:
  - The grid materials match Python **bit for bit**: copper fractions and outline mask exactly, conductivities
    and via areas to 1e-14.
  - Rust Jacobi and Rust AMG both agree with Python to 1e-6 of each field's scale (the measured worst
    difference is about 1e-9).
  - The report text is identical.
  - Jacobi iteration counts equal Python's, and AMG takes at most twice pyamg's.
- `tests/analytic.rs` checks plate convection against P/(2hA) and track IR drop against ρL/(wt).
- `tests/determinism.rs` checks that results are bit-identical on 1 and 4 threads.

A 1 % change to FR4's through-plane conductivity, or dropping convection from the bottom face, fails every
regression case.

**Real board.** A 6-layer, 170 × 215 mm power-distribution board went through the full workflow of a pre-fab
power-integrity check: 23 power paths solved at 0.1 mm (1701 × 2151 × 13 = 47.6 M cells), 4 ground paths and one
whole-board thermal solve at 0.25 mm (7.6 M cells), 24 heat sources plus the Joule heat of every path. Both tools
ran on the same exported geometry. Compared with the Python original:

- All 27 paths gave the same pass/fail result and the same location of peak current density. All 24 component
  temperatures, the board maximum and the hottest cell were identical.
- The 13 paths that `auto` solves with Jacobi matched Python exactly. The 14 that use AMG differed by at most
  one reported digit (0.1 mV on a 62 mV drop, 0.3 A/mm² on 636 A/mm²).
- Converged to a 1e-13 residual, Jacobi and AMG both land on the AMG value. At the default 1e-9 tolerance it is
  Python's Jacobi-CG that has not fully converged on a grid this size; AMG has.

## Performance

**Real board** (above), Intel i5-1240P, 16 threads:

| Step | Python original | kicad-thermal (`auto`) |
|---|---|---|
| Full power-integrity check (27 DC paths + thermal) | 63 min | 128 s |
| Large DC path at 0.1 mm (AMG) | 255–393 s | 7.5 s |
| Ground path at 0.25 mm (AMG) | 290–330 s | 7–8 s |
| Thermal solve at 0.25 mm, 7.6 M cells | 734 s | 26 s |

The Python run shared the machine with other work, so treat the ratios (28–44 times per step) as approximate.
Peak memory for the full kicad-thermal run was 3.7 GB.

**Synthetic boards**, thermal solve, one CPU core (seconds):

| Board | Cells | Python Jacobi | Python + pyamg | Rust Jacobi | Rust AMG |
|---|---|---|---|---|---|
| 2-layer 60×40 mm, 0.25 mm | 0.23M | 4.4 | 3.4 | 3.4 | 1.2 |
| 2-layer strip 200×30 mm, 0.25 mm | 0.58M | 23.0 | 10.0 | 16.6 | 4.7 |
| 4-layer 120×90 mm, 0.25 mm | 1.56M | 88.0 | 24.8 | 81.8 | 12.5 |
| 6-layer heavy copper 100×80 mm, 0.25 mm | 1.66M | 88.9 | 29.3 | 82.7 | 15.0 |
| 2-layer 60×40 mm, 0.1 mm | 1.44M | 96.1 | 34.1 | 93.2 | 16.8 |

Most of the gain comes from the algorithm, not the language. The Jacobi solver is limited by memory bandwidth,
so on one core the Rust version is only slightly faster than SciPy's. AMG cuts CG iterations from thousands to
about 20–90. Vector operations and SpMV run in parallel with rayon, but the AMG smoother (symmetric Gauss–Seidel)
is sequential, so on the real board the AMG path keeps only about 3 cores busy.

## Known differences from kicad_fdtd

- Copper sampling uses scanline intervals instead of shapely's exact predicates. A sample point lying exactly
  on a polygon edge could be classified differently.
- A polygon is treated as its outline minus the union of its holes, which matches shapely's `buffer(0)` repair
  of overlapping holes. A self-intersecting outline is filled even–odd, which `buffer(0)` may repair differently.
- AMG follows pyamg's smoothed-aggregation defaults (θ = 0, standard aggregation, ω = 4/3, symmetric
  Gauss–Seidel, Galerkin coarse operators) but differs in spectral-radius estimation and coarse-grid size. It
  converges to the same answer, not through the same iterates.

## Development

```bash
cargo test --release --locked          # the regression suite is the contract
```

The Python original is not part of this repository. `scripts/fetch_kicad_fdtd.sh` clones it at `74e21eb` into
`third_party/kicad_fdtd` (git-ignored) and applies `patches/kicad_fdtd-amg-preconditioner.patch`, which adds
pyamg's AMG and a `precond` option so both solvers can be compared. The scripts below default to that location;
set `KICAD_FDTD` to use another checkout. They need numpy, scipy and shapely, plus pyamg for AMG.

| Task | Command |
|---|---|
| Compare both tools on any geometry JSON (wall times, worst field differences, report text) | `python scripts/compare_live.py board_geometry.json -- --power U3=1.8 --h 5 --tamb 45 --net /12V --from U3.5 --to J7.1 --current 5` |
| Benchmark Rust against Python on synthetic boards, at several thread counts | `python bench/boards.py && python bench/rust_vs_python.py` |
| Diff grid materials cell by cell (debugging a model change) | `cargo build --release --examples && python bench/cmp_grid.py B_2L_compact 0.25` |
| Regenerate the Python reference (only when the Python original changes) | `python scripts/gen_reference.py` (and `scripts/make_fixtures.py` if the fixtures change) |

`compare_live.py` exits non-zero if any normalized difference exceeds `--tol` (default 1e-6).

| Path | Contents |
|---|---|
| `src/geometry.rs` | Exporter JSON, stackup, model frame |
| `src/raster.rs` | Point-in-polygon-union sampling |
| `src/grid.rs` | Grid, copper fractions, via barrels, conductivities |
| `src/solve.rs` | Conductance links, DC solve, thermal solve |
| `src/cg.rs` | PCG, step for step as `scipy.sparse.linalg.cg` |
| `src/amg.rs` | Smoothed-aggregation AMG |
| `src/sparse.rs` | CSR matrix and deterministic parallel kernels |
| `src/report.rs` | Report text and result JSON |
| `examples/dump_grid.rs` | Writes the grid materials as raw f64, for `bench/cmp_grid.py` |
| `tests/` | Regression, analytic and determinism tests; fixtures and Python reference data |
| `scripts/` | The geometry exporter (`export_geometry.py`); fetching and patching kicad_fdtd, reference generation, live comparison, the KiCad 9 shim |
| `bench/` | Synthetic benchmark boards and the Rust vs Python benchmark |
| `patches/` | The AMG patch for kicad_fdtd |

## Roadmap

- A batch mode: one grid per cell size, many DC paths, and one thermal solve with their summed Joule heat, from
  the command line (the library already supports this). Reusing the AMG hierarchy across solves on the same
  matrix would cut setup time further.
- A parallel AMG smoother (l1-Jacobi or Chebyshev), so the AMG path scales with cores.
- A matrix-free stencil kernel for the regular grid, which roughly halves memory traffic per iteration.

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT License ([LICENSE-MIT](LICENSE-MIT))

at your option.

kicad-thermal is derived from [kicad_fdtd](https://github.com/Jerome-Graves/kicad_fdtd), Copyright (c) 2026
Jerome Graves, under the MIT License. Its copyright and permission notice are reproduced in [NOTICE](NOTICE) and
must be kept under either license.

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in this work by
you, as defined in the Apache-2.0 license, shall be dual licensed as above, without any additional terms or
conditions.
