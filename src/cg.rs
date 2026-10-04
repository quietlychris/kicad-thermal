//! Preconditioned conjugate gradient, step for step as `scipy.sparse.linalg.cg`:
//! zero start, stop when ||r|| < tol * ||b|| (checked at the top of each iteration on the
//! true residual), `maxiter` iterations at most.

use crate::amg::{Amg, AmgInfo, AmgOptions};
use crate::sparse::{chunked_sum, norm, Csr, CHUNK};
use anyhow::{bail, Result};
use rayon::prelude::*;
use std::time::Instant;

/// Unknowns below which `Precond::Auto` uses Jacobi (AMG setup costs more than it saves).
pub const AMG_MIN_CELLS: usize = 50_000;

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum Precond {
    /// AMG for systems of at least AMG_MIN_CELLS unknowns, Jacobi below
    Auto,
    /// smoothed-aggregation algebraic multigrid
    Amg,
    /// diagonal (the original kicad_fdtd solver)
    Jacobi,
}

#[derive(Clone, Debug)]
pub struct SolveStats {
    pub precond: &'static str,
    pub n: usize,
    pub iterations: usize,
    pub setup_s: f64,
    pub solve_s: f64,
    /// AMG hierarchy sizes, when AMG was used
    pub amg: Option<AmgInfo>,
    /// set when AMG setup failed and Jacobi was used instead
    pub fallback: Option<String>,
}

pub const DEFAULT_TOL: f64 = 1e-9;
pub const DEFAULT_MAXITER: usize = 20_000;

/// Solve A x = b for SPD A.
pub fn solve(
    a: &Csr,
    b: &[f64],
    tol: f64,
    maxiter: usize,
    precond: Precond,
) -> Result<(Vec<f64>, SolveStats)> {
    let n = a.nrows;
    let want_amg = precond == Precond::Amg || (precond == Precond::Auto && n >= AMG_MIN_CELLS);
    let t0 = Instant::now();
    let mut fallback = None;
    let amg = if want_amg {
        match Amg::new(a, &AmgOptions::default()) {
            Ok(m) => Some(m),
            Err(e) if precond == Precond::Auto => {
                fallback = Some(format!("AMG setup failed ({e}); used Jacobi"));
                None
            }
            Err(e) => return Err(e),
        }
    } else {
        None
    };
    let dinv: Vec<f64> = if amg.is_none() {
        a.diagonal().iter().map(|d| 1.0 / d).collect()
    } else {
        Vec::new()
    };
    let setup_s = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let (x, iterations, converged) = pcg(a, b, tol, maxiter, |r, z| match &amg {
        Some(m) => {
            m.apply(r, z);
            crate::sparse::dot(r, z)
        }
        None => jacobi_apply_dot(&dinv, r, z),
    });
    let stats = SolveStats {
        precond: if amg.is_some() { "amg" } else { "jacobi" },
        n,
        iterations,
        setup_s,
        solve_s: t1.elapsed().as_secs_f64(),
        amg: amg.as_ref().map(|m| m.info.clone()),
        fallback,
    };
    if !converged {
        bail!("CG did not converge (info {})", maxiter);
    }
    Ok((x, stats))
}

/// z = dinv * r, returning r . z
fn jacobi_apply_dot(dinv: &[f64], r: &[f64], z: &mut [f64]) -> f64 {
    let parts: Vec<f64> = z
        .par_chunks_mut(CHUNK)
        .enumerate()
        .map(|(c, zs)| {
            let r0 = c * CHUNK;
            let mut acc = 0.0;
            for (o, zi) in zs.iter_mut().enumerate() {
                let i = r0 + o;
                *zi = r[i] * dinv[i];
                acc += r[i] * *zi;
            }
            acc
        })
        .collect();
    parts.iter().sum()
}

/// PCG with a preconditioner callback `m(r, z) -> r . z`.
/// Returns (x, iterations, converged).
pub fn pcg(
    a: &Csr,
    b: &[f64],
    tol: f64,
    maxiter: usize,
    m: impl Fn(&[f64], &mut [f64]) -> f64,
) -> (Vec<f64>, usize, bool) {
    let n = b.len();
    let bnrm2 = norm(b);
    let mut x = vec![0.0; n];
    if bnrm2 == 0.0 {
        return (x, 0, true);
    }
    let atol = tol * bnrm2;
    let mut r = b.to_vec();
    let mut z = vec![0.0; n];
    let mut p = vec![0.0; n];
    let mut q = vec![0.0; n];
    let mut rr = chunked_sum(n, |rg| rg.map(|i| r[i] * r[i]).sum());
    let mut rho_prev = 0.0;
    for it in 0..maxiter {
        if rr.sqrt() < atol {
            return (x, it, true);
        }
        let rho = m(&r, &mut z);
        if it > 0 {
            let beta = rho / rho_prev;
            p.par_iter_mut()
                .zip(z.par_iter())
                .for_each(|(pi, zi)| *pi = *pi * beta + *zi);
        } else {
            p.copy_from_slice(&z);
        }
        let pq = a.spmv_dot(&p, &mut q);
        let alpha = rho / pq;
        let parts: Vec<f64> = x
            .par_chunks_mut(CHUNK)
            .zip(r.par_chunks_mut(CHUNK))
            .enumerate()
            .map(|(c, (xs, rs))| {
                let r0 = c * CHUNK;
                let mut acc = 0.0;
                for o in 0..xs.len() {
                    let i = r0 + o;
                    xs[o] += alpha * p[i];
                    rs[o] -= alpha * q[i];
                    acc += rs[o] * rs[o];
                }
                acc
            })
            .collect();
        rr = parts.iter().sum();
        rho_prev = rho;
    }
    (x, maxiter, false)
}
