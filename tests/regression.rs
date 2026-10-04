//! Accuracy regression suite against the Python reference (kicad_fdtd + AMG patch).
//!
//! For every case in tests/fixtures/cases.json:
//! 1. the grid materials must match Python exactly: copper fractions and the outline mask
//!    bit for bit, conductivities and via areas to 1e-14 relative;
//! 2. the Rust solve with BOTH preconditioners must match BOTH Python solves (the original
//!    Jacobi solver and the AMG one) to `FIELD_TOL` of the field's scale: temperature to
//!    1e-6 of the maximum rise, voltage to 1e-6 of the drop, current density to 1e-6 of
//!    its maximum. CG stops at a relative residual of 1e-9, so two correct solves agree to
//!    ~1e-9; a real model change (materials, links, sources) moves results by far more;
//! 3. the report text must match Python's exactly;
//! 4. iteration counts guard solver quality: Rust Jacobi within 10% of Python's (same
//!    algorithm, different summation order), Rust AMG at most twice pyamg's.

mod common;
use common::*;
use kicad_thermal::{Grid, Precond};

const FIELD_TOL: f64 = 1e-6;
const SCALAR_TOL: f64 = 1e-6;

fn check_materials(grid: &Grid, r: &Reference) {
    let (nx, ny) = (grid.nx, grid.ny);
    let g = &r.meta["grid"];
    assert_eq!(
        (nx, ny, grid.nz),
        (
            g["nx"].as_u64().unwrap() as usize,
            g["ny"].as_u64().unwrap() as usize,
            g["nz"].as_u64().unwrap() as usize
        ),
        "grid shape"
    );
    for (l, f) in grid.cu_frac.iter().enumerate() {
        let py = r.grid_array(&format!("cu_frac_{l}"));
        let bad = f.iter().zip(&py).filter(|(a, b)| a != b).count();
        assert_eq!(
            bad, 0,
            "copper fraction differs from Python in {bad} cells on layer {}",
            grid.bd.layers[l]
        );
    }
    let inside = r.grid_array("inside");
    let bad = grid
        .inside
        .iter()
        .zip(&inside)
        .filter(|(a, b)| (**a as u8 as f64) != **b)
        .count();
    assert_eq!(bad, 0, "outline mask differs in {bad} cells");
    for (key, ours) in [
        ("kx", grid.kx.clone()),
        ("kz", grid.kz.clone()),
        ("via_area", grid.via_copper_area(None)),
    ] {
        let py = r.grid_array(key);
        let worst = ours
            .iter()
            .zip(&py)
            .map(|(a, b)| rel(*a, *b))
            .fold(0.0, f64::max);
        assert!(
            worst <= 1e-14,
            "{key} differs from Python by {worst:e} relative"
        );
    }
}

/// Returns the worst normalized field difference found.
fn check_against(run: &Run, r: &Reference, py: &str, label: &str) -> f64 {
    let res = r.result(py);
    let case = r.case();
    let tamb = case["tamb"].as_f64().unwrap();
    // thermal
    let t_ref = r.result_array(py, "T").unwrap();
    let py_tmax = res["thermal"]["Tmax"].as_f64().unwrap();
    let rise = (py_tmax - tamb).abs().max(1e-12);
    let dt = max_abs_diff(&run.th.t, &t_ref);
    assert!(
        dt <= FIELD_TOL * rise,
        "{label} vs python-{py}: max |dT| = {dt:e} C = {:e} of the {rise:.3} K rise",
        dt / rise
    );
    let mut worst = dt / rise;
    assert!(
        (run.th.tmax - py_tmax).abs() <= FIELD_TOL * rise,
        "{label} vs python-{py}: Tmax {} vs {}",
        run.th.tmax,
        py_tmax
    );
    let pt = res["thermal"]["P_total"].as_f64().unwrap();
    assert!(
        rel(run.th.p_total, pt) <= 1e-12 || (run.th.p_total - pt).abs() < 1e-15,
        "{label}: P_total {} vs {}",
        run.th.p_total,
        pt
    );
    // DC
    match (&run.dc, res.get("dc")) {
        (Some(dc), Some(rd)) => {
            assert_eq!(
                dc.cells as u64,
                rd["cells"].as_u64().unwrap(),
                "{label}: DC cell count"
            );
            for (name, ours) in [
                ("drop_V", dc.drop_v),
                ("R_ohm", dc.r_ohm),
                ("joule_W", dc.joule_w),
                ("Jmax", dc.jmax),
            ] {
                let theirs = rd[name].as_f64().unwrap();
                assert!(
                    rel(ours, theirs) <= SCALAR_TOL,
                    "{label} vs python-{py}: {name} {ours:e} vs {theirs:e}"
                );
            }
            let v_ref = r.result_array(py, "V").unwrap();
            let mask_bad =
                dc.v.iter()
                    .zip(&v_ref)
                    .filter(|(a, b)| a.is_nan() != b.is_nan())
                    .count();
            assert_eq!(
                mask_bad, 0,
                "{label}: DC copper mask differs in {mask_bad} cells"
            );
            let dv =
                dc.v.iter()
                    .zip(&v_ref)
                    .filter(|(a, _)| !a.is_nan())
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max);
            let drop = rd["drop_V"].as_f64().unwrap().abs();
            assert!(
                dv <= FIELD_TOL * drop,
                "{label} vs python-{py}: max |dV| = {:e} of the drop",
                dv / drop
            );
            let j_ref = r.result_array(py, "J").unwrap();
            let dj = max_abs_diff(&dc.j, &j_ref);
            let jmax = rd["Jmax"].as_f64().unwrap();
            assert!(
                dj <= FIELD_TOL * jmax,
                "{label} vs python-{py}: max |dJ| = {:e} of Jmax",
                dj / jmax
            );
            worst = worst.max(dv / drop).max(dj / jmax);
        }
        (None, None) => {}
        _ => panic!("{label}: DC presence differs from the reference"),
    }
    worst
}

fn check_case(name: &str) {
    let (bd, r) = load_case(name);
    let case = r.case().clone();
    let grid = Grid::new(&bd, case["cell"].as_f64().unwrap()).unwrap();
    check_materials(&grid, &r);
    for (label, pc) in [("rust-jacobi", Precond::Jacobi), ("rust-amg", Precond::Amg)] {
        let run = run_case(&grid, &case, pc);
        let worst =
            check_against(&run, &r, "jacobi", label).max(check_against(&run, &r, "amg", label));
        assert_eq!(
            run.report,
            r.result("jacobi")["report"].as_str().unwrap(),
            "{label}: report text differs from Python"
        );
        // iteration counts
        let py_it = |p: &str, what: &str| r.result(p)[what]["iterations"].as_u64().unwrap() as f64;
        let ours = run.th.stats.iterations as f64;
        match pc {
            Precond::Jacobi => {
                let theirs = py_it("jacobi", "thermal");
                assert!(
                    (ours - theirs).abs() <= 0.10 * theirs,
                    "{label}: {ours} thermal iterations vs Python's {theirs}"
                );
            }
            _ => {
                let theirs = py_it("amg", "thermal");
                assert!(
                    ours <= 2.0 * theirs,
                    "{label}: {ours} thermal iterations vs pyamg's {theirs}"
                );
            }
        }
        eprintln!(
            "{name:<15} {label:<11} worst field diff {worst:.1e} | thermal {:>4} it ({}){}",
            run.th.stats.iterations,
            run.th.stats.precond,
            run.dc.as_ref().map_or(String::new(), |d| format!(
                ", DC {:>4} it ({})",
                d.stats.iterations, d.stats.precond
            ))
        );
    }
}

macro_rules! cases {
    ($($name:ident),* $(,)?) => { $( #[test] fn $name() { check_case(stringify!($name)); } )* };
}

cases!(
    a4_power,
    b2_compact,
    c6_heavy,
    d2_strip,
    e2_lshape_vias,
    plate,
    track
);

/// Every case in cases.json has a test above (catches cases added without one).
#[test]
fn all_cases_covered() {
    let cases: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(root().join("fixtures/cases.json")).unwrap())
            .unwrap();
    let names: Vec<&str> = cases
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c["name"].as_str().unwrap())
        .collect();
    assert_eq!(
        names,
        [
            "a4_power",
            "b2_compact",
            "c6_heavy",
            "d2_strip",
            "e2_lshape_vias",
            "plate",
            "track"
        ]
    );
}
