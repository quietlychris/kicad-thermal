//! Steady-state thermal and DC IR-drop solver for KiCad boards.
//!
//! A Rust port of the thermal module of [kicad_fdtd](https://github.com/Jerome-Graves/kicad_fdtd)
//! (MIT, Jerome Graves): the same finite-volume model on the board's exact copper, solved
//! with conjugate gradient and either the original Jacobi preconditioner or smoothed-
//! aggregation algebraic multigrid.
//!
//! Pipeline: [`geometry::Board::load`] (the exporter's JSON) → [`grid::Grid::new`] →
//! [`solve::solve_dc`] (optional) → [`solve::solve_thermal`] → [`report`].

pub mod amg;
pub mod cg;
pub mod geometry;
pub mod grid;
pub mod raster;
pub mod report;
pub mod solve;
pub mod sparse;

pub use cg::Precond;
pub use geometry::Board;
pub use grid::Grid;
pub use solve::{solve_dc, solve_thermal, DcResult, ThermalResult};
