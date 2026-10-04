# kicad-thermal — project brief for Claude Code

`kicad-thermal` is a Rust port of the steady-state thermal and DC IR-drop solver in
[kicad_fdtd](https://github.com/Jerome-Graves/kicad_fdtd) (Python, MIT). The work is the port
itself and proving it matches the original: the same model, results agreeing with Python to
solver tolerance, and faster through an AMG preconditioner. The Python original is not in this repo:
`scripts/fetch_kicad_fdtd.sh` fetches it (pinned commit + the AMG patch) into `third_party/kicad_fdtd` when it is
needed, as the test reference and as the exporter of a `.kicad_pcb` to the geometry JSON the Rust tool reads (that
export needs KiCad's `pcbnew`).

Ask the maintainer before changing anything listed under "Ask first".

## Toolchain

Rust `1.99.0`. Every build and test uses `--locked`; `Cargo.lock` is committed and only changes on purpose.
The Python reference (only for regenerating references, comparisons and benchmarks) runs on numpy, scipy, shapely
and pyamg; the versions used for the committed reference are recorded in `tests/reference/` metadata.

## Where things are

| Path | What |
|---|---|
| `README.md` | Usage, validation (tests and the real-board comparison), performance, known differences, layout |
| `src/` | The solver; the module map is in `README.md` |
| `tests/regression.rs` | 7 fixture boards against the Python reference data in `tests/reference/` |
| `tests/analytic.rs`, `tests/determinism.rs` | Physics checks, and bit-identical results on 1 vs 4 threads |
| `scripts/fetch_kicad_fdtd.sh` | Clones kicad_fdtd at `74e21eb` into `third_party/` and applies `patches/` |
| `scripts/compare_live.py` | Runs Python and Rust on any geometry JSON and diffs the results |
| `scripts/gen_reference.py`, `scripts/make_fixtures.py` | Regenerate the Python reference and fixtures (see the rules below) |
| `scripts/kicad9_shim/` | Runtime shim so kicad_fdtd's exporter works on KiCad 9 |
| `bench/` | Synthetic benchmark boards, Rust vs Python benchmark, grid-materials diff |
| `NOTICE`, `LICENSE-MIT`, `LICENSE-APACHE` | Dual license; kicad_fdtd's MIT notice must stay in `NOTICE` |

## Status

- **The port is complete for thermal and DC.** All 16 tests pass; grid materials are bit-identical to Python,
  fields agree to ~1e-9 against a 1e-6 tolerance, report text is identical.
- **Real-board comparison done** (6-layer 170 × 215 mm board, 47.6 M cells at 0.1 mm, 27 DC paths + thermal, via the
  library): identical pass/fail, peak locations and temperatures; Jacobi paths exact; AMG paths within one reported
  digit and converged further than Python's Jacobi at 1e-9. About 30× faster (128 s vs 63 min). See `README.md`.

## Next steps

1. **Batch mode** (CLI): one grid per cell size, many DC paths, one thermal solve with their summed Joule heat, from a
   JSON case file. Add a test that batched and individual results are bit-identical. Reuse the AMG hierarchy across
   solves on the same matrix.
2. **Parallel AMG smoother** (l1-Jacobi or Chebyshev), keeping the V-cycle symmetric. The "AMG iterations <= 2x pyamg"
   check must still pass. Benchmark at 1 and N threads.
3. **Matrix-free stencil kernel** for the Jacobi path, which should roughly halve memory traffic per iteration.
4. **Optional**: clear the clippy warnings without changing results.

## Before every commit

```bash
cargo fmt --all
cargo clippy --all-targets --locked
cargo test --locked --release
```

## Working rules

- **The regression suite is the contract.** Never loosen `FIELD_TOL`/`SCALAR_TOL`, the bit-exact
  materials check, or the iteration-count checks to make a change pass.
- **Never regenerate `tests/reference/` to fix a Rust failure.** References come only from the
  Python original, only when the Python itself changes, and only with the maintainer's OK.
- **After any model change**, inject a small deliberate error once to confirm the suite still
  catches it, then revert.
- **Keep results deterministic.** Parallel reductions use fixed `CHUNK`s combined in order. Don't
  use rayon's unordered reductions in solver paths.
- **The Python baseline is `kicad_fdtd` with `precond="jacobi"`**, which is exactly the original
  algorithm. Every comparison or benchmark uses the same board, cell size, powers and current for
  both tools, and states the core count.
- **Copper sampling must match shapely.** A polygon is its outline (even-odd) minus the union of its
  holes, and samples count only when strictly inside. Even-odd across holes was a real bug, caught
  on overlapping plane cutouts.
- **The model's limits are kicad_fdtd's**, and the port must not "fix" them silently: there is no
  enclosure or radiation, and peak temperature and neck voltage drop depend on cell size. Any
  deliberate model change is a divergence from the reference (see Ask first).
- **Keeping the patch file current**: changes to the Python original are made in `third_party/kicad_fdtd` and kept
  as a diff against `74e21eb`. Refresh `patches/` with `git -C third_party/kicad_fdtd diff` after any change.
- **Licensing**: dual MIT OR Apache-2.0; `NOTICE` and `LICENSE-MIT` must keep Jerome Graves' MIT notice.

## Ask first

- Any change to the physical model or constants, or any deliberate divergence from kicad_fdtd.
- Loosening a tolerance, or regenerating the reference data.
- New dependencies, an edition or MSRV bump, or `unsafe`.
- Publishing: crates.io, a public repo, or an upstream pull request with the AMG patch.
