//! Physics checks independent of the Python reference (ported from kicad_fdtd's
//! tests/test_thermal.py), run with both preconditioners.

mod common;
use common::*;
use kicad_thermal::{grid::RHO_CU, solve_dc, solve_thermal, Board, Grid, Precond};

/// A copper-clad plate with P watts and convection h on both faces settles at
/// T_amb + P / (2 h A): mean rise within 3 %, spread small.
#[test]
fn plate_convection() {
    let bd = Board::load(fixture("plate")).unwrap();
    let g = Grid::new(&bd, 0.5).unwrap();
    let (p, h) = (1.0, 10.0);
    let expect = p / (2.0 * h * 20e-3 * 8e-3);
    for pc in [Precond::Jacobi, Precond::Amg] {
        let th = solve_thermal(&g, &[("P".into(), p)], h, 25.0, None, pc).unwrap();
        let k = g.k_cu[0];
        let top: Vec<f64> = (0..g.nx * g.ny).map(|c| th.t[c * g.nz + k]).collect();
        let rise = top.iter().sum::<f64>() / top.len() as f64 - 25.0;
        let spread = th.t.iter().cloned().fold(f64::MIN, f64::max)
            - th.t.iter().cloned().fold(f64::MAX, f64::min);
        assert!(
            (rise - expect).abs() / expect < 0.03,
            "{pc:?}: mean rise {rise} vs {expect}"
        );
        assert!(spread < 0.15 * expect, "{pc:?}: spread {spread}");
    }
}

/// A 0.2 mm x 35 um track carrying 1 A over 11.4 mm drops I rho L / (w t) = 27.4 mV;
/// Joule heat equals I^2 R.
#[test]
fn track_ir_drop() {
    let bd = Board::load(fixture("track")).unwrap();
    let g = Grid::new(&bd, 0.1).unwrap();
    let (l, w, t) = (11.4e-3, 0.2e-3, 35e-6);
    let expect = RHO_CU * l / (w * t);
    for pc in [Precond::Jacobi, Precond::Amg] {
        let dc = solve_dc(&g, "/T", ("P", "1"), ("P", "2"), 1.0, pc).unwrap();
        assert!(
            (dc.drop_v - expect).abs() / expect < 0.1,
            "{pc:?}: drop {} vs {}",
            dc.drop_v,
            expect
        );
        let jexp = 1.0 / (w * t) / 1e6;
        assert!(
            (dc.jmax - jexp).abs() / jexp < 0.25,
            "{pc:?}: Jmax {} vs {}",
            dc.jmax,
            jexp
        );
        assert!(
            (dc.joule_w - dc.drop_v).abs() / dc.drop_v < 0.02,
            "{pc:?}: Joule {} vs I*V {}",
            dc.joule_w,
            dc.drop_v
        );
    }
}

/// Errors that kicad_fdtd raises must surface as errors, not panics or wrong answers.
#[test]
fn bad_inputs_are_errors() {
    let bd = Board::load(fixture("track")).unwrap();
    let g = Grid::new(&bd, 0.1).unwrap();
    assert!(solve_thermal(
        &g,
        &[("NOPE".into(), 1.0)],
        10.0,
        25.0,
        None,
        Precond::Jacobi
    )
    .is_err());
    // pad on a net with no copper under it
    assert!(solve_dc(&g, "/CLK_P", ("P", "1"), ("P", "2"), 1.0, Precond::Jacobi).is_err());
}
