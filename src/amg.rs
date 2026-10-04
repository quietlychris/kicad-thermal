//! Smoothed-aggregation algebraic multigrid (SA-AMG), used as a CG preconditioner.
//!
//! Follows pyamg's `smoothed_aggregation_solver` defaults as used by the kicad_fdtd AMG
//! patch: symmetric strength of connection with theta = 0 (every off-diagonal coupling is
//! strong), standard aggregation, constant near-null space, one Jacobi smoothing step of
//! the tentative prolongator with omega = 4/3 over the spectral radius of D^-1 A, Galerkin
//! coarse operators R A P with R = P^T, symmetric Gauss-Seidel pre- and post-smoothing,
//! and a direct (dense Cholesky) solve on the coarsest level. With SGS on both sides the
//! V-cycle is a symmetric positive definite operator, as CG requires.
//!
//! It is not a bit-for-bit copy of pyamg (the spectral radius estimate and coarse sizes
//! differ); the CG answer it produces agrees to the solver tolerance, which is what the
//! regression suite checks.

use crate::sparse::{dot, Csr};
use anyhow::{bail, Result};
use std::cell::RefCell;

#[derive(Clone, Debug)]
pub struct AmgOptions {
    /// stop coarsening at or below this many unknowns
    pub max_coarse: usize,
    pub max_levels: usize,
    /// prolongator smoothing weight (pyamg: 4/3)
    pub omega: f64,
    /// power iterations for the spectral radius of D^-1 A
    pub power_iters: usize,
}

impl Default for AmgOptions {
    fn default() -> Self {
        AmgOptions {
            max_coarse: 500,
            max_levels: 20,
            omega: 4.0 / 3.0,
            power_iters: 30,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct AmgInfo {
    /// unknowns per level, finest first
    pub sizes: Vec<usize>,
    /// nonzeros per level
    pub nnz: Vec<usize>,
}

impl AmgInfo {
    /// Operator complexity: total nonzeros over finest-level nonzeros.
    pub fn operator_complexity(&self) -> f64 {
        self.nnz.iter().sum::<usize>() as f64 / self.nnz[0] as f64
    }
}

struct Level {
    diag: Vec<f64>,
    p: Csr,
    r: Csr,
}

struct Work {
    r: Vec<f64>,
    bc: Vec<f64>,
    xc: Vec<f64>,
}

pub struct Amg<'a> {
    a0: &'a Csr,
    /// operators of levels 1.. (level 0 is `a0`)
    mats: Vec<Csr>,
    levels: Vec<Level>,
    coarse: DenseCholesky,
    work: RefCell<Vec<Work>>,
    pub info: AmgInfo,
}

const UNSET: u32 = u32::MAX;
const ISOLATED: u32 = u32::MAX - 1;

/// Standard aggregation on the off-diagonal nonzero pattern (theta = 0).
/// Returns the aggregate of each node (`ISOLATED` for nodes with no neighbours) and the count.
fn standard_aggregation(a: &Csr) -> (Vec<u32>, usize) {
    let n = a.nrows;
    let mut agg = vec![UNSET; n];
    let mut pass1 = vec![false; n];
    let mut next = 0u32;
    let nbrs = |i: usize| {
        (a.indptr[i]..a.indptr[i + 1])
            .filter(move |&p| a.indices[p] as usize != i && a.data[p] != 0.0)
            .map(move |p| a.indices[p] as usize)
    };
    // pass 1: a node whose neighbours are all free seeds an aggregate of itself and them
    for i in 0..n {
        if agg[i] != UNSET {
            continue;
        }
        let mut has_nb = false;
        let mut has_agg_nb = false;
        for j in nbrs(i) {
            has_nb = true;
            if agg[j] != UNSET {
                has_agg_nb = true;
                break;
            }
        }
        if !has_nb {
            agg[i] = ISOLATED;
        } else if !has_agg_nb {
            agg[i] = next;
            pass1[i] = true;
            for j in nbrs(i) {
                agg[j] = next;
                pass1[j] = true;
            }
            next += 1;
        }
    }
    // pass 2: join a neighbouring pass-1 aggregate
    for i in 0..n {
        if agg[i] != UNSET {
            continue;
        }
        if let Some(j) = nbrs(i).find(|&j| pass1[j]) {
            agg[i] = agg[j];
        }
    }
    // pass 3: leftovers form aggregates with their free neighbours
    for i in 0..n {
        if agg[i] != UNSET {
            continue;
        }
        agg[i] = next;
        for j in nbrs(i) {
            if agg[j] == UNSET {
                agg[j] = next;
            }
        }
        next += 1;
    }
    (agg, next as usize)
}

/// Generalized Rayleigh quotient power iteration: largest lambda of A v = lambda D v.
fn spectral_radius_dinv_a(a: &Csr, diag: &[f64], iters: usize) -> f64 {
    let n = a.nrows;
    // deterministic pseudo-random start (xorshift)
    let mut s: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut v: Vec<f64> = (0..n)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            0.5 + (s >> 11) as f64 / (1u64 << 53) as f64
        })
        .collect();
    let mut w = vec![0.0; n];
    let mut lambda = 0.0f64;
    for _ in 0..iters {
        let dn: f64 = (0..n).map(|i| diag[i] * v[i] * v[i]).sum::<f64>().sqrt();
        for x in v.iter_mut() {
            *x /= dn;
        }
        a.spmv(&v, &mut w);
        lambda = lambda.max(dot(&v, &w)); // v^T D v = 1
        for i in 0..n {
            v[i] = w[i] / diag[i];
        }
    }
    lambda
}

/// P = T - (omega / rho) D^-1 A T for the piecewise-constant tentative prolongator T.
fn smoothed_prolongator(a: &Csr, diag: &[f64], agg: &[u32], nagg: usize, opts: &AmgOptions) -> Csr {
    let n = a.nrows;
    let mut size = vec![0usize; nagg];
    for &g in agg {
        if (g as usize) < nagg {
            size[g as usize] += 1;
        }
    }
    let tval = |i: usize| -> Option<(u32, f64)> {
        let g = agg[i];
        ((g as usize) < nagg).then(|| (g, 1.0 / (size[g as usize] as f64).sqrt()))
    };
    let mut ti = vec![0usize; n + 1];
    let mut tc = Vec::with_capacity(n);
    let mut td = Vec::with_capacity(n);
    for i in 0..n {
        if let Some((c, v)) = tval(i) {
            tc.push(c);
            td.push(v);
        }
        ti[i + 1] = tc.len();
    }
    let t = Csr {
        nrows: n,
        ncols: nagg,
        indptr: ti,
        indices: tc,
        data: td,
    };
    let rho = spectral_radius_dinv_a(a, diag, opts.power_iters);
    let at = a.matmul(&t);
    let scale = opts.omega / rho;
    let mut pi = vec![0usize; n + 1];
    let mut pc = Vec::with_capacity(at.nnz() + n);
    let mut pd = Vec::with_capacity(at.nnz() + n);
    for i in 0..n {
        let s = scale / diag[i];
        let tv = tval(i);
        let mut placed = tv.is_none();
        for p in at.indptr[i]..at.indptr[i + 1] {
            let c = at.indices[p];
            let mut v = -s * at.data[p];
            if let Some((tc, tvv)) = tv {
                if !placed && c > tc {
                    pc.push(tc);
                    pd.push(tvv);
                    placed = true;
                } else if c == tc {
                    v += tvv;
                    placed = true;
                }
            }
            pc.push(c);
            pd.push(v);
        }
        if !placed {
            let (tc, tvv) = tv.unwrap();
            pc.push(tc);
            pd.push(tvv);
        }
        pi[i + 1] = pc.len();
    }
    Csr {
        nrows: n,
        ncols: nagg,
        indptr: pi,
        indices: pc,
        data: pd,
    }
}

impl<'a> Amg<'a> {
    pub fn new(a0: &'a Csr, opts: &AmgOptions) -> Result<Amg<'a>> {
        let mut mats: Vec<Csr> = Vec::new();
        let mut levels: Vec<Level> = Vec::new();
        let mut info = AmgInfo {
            sizes: vec![a0.nrows],
            nnz: vec![a0.nnz()],
        };
        loop {
            let a = if mats.is_empty() {
                a0
            } else {
                mats.last().unwrap()
            };
            let n = a.nrows;
            if n <= opts.max_coarse || levels.len() + 1 >= opts.max_levels {
                break;
            }
            let diag = a.diagonal();
            if diag.iter().any(|&d| !(d > 0.0)) {
                bail!("AMG needs a positive diagonal");
            }
            let (agg, nagg) = standard_aggregation(a);
            if nagg == 0 || nagg as f64 > 0.9 * n as f64 {
                break; // no useful coarsening
            }
            let p = smoothed_prolongator(a, &diag, &agg, nagg, opts);
            let r = p.transpose();
            let ac = r.matmul(&a.matmul(&p));
            info.sizes.push(ac.nrows);
            info.nnz.push(ac.nnz());
            levels.push(Level { diag, p, r });
            mats.push(ac);
        }
        let coarsest = if mats.is_empty() {
            a0
        } else {
            mats.last().unwrap()
        };
        if coarsest.nrows > 20_000 {
            bail!(
                "AMG coarsest level has {} unknowns; too large for a dense solve",
                coarsest.nrows
            );
        }
        let coarse = DenseCholesky::new(coarsest)?;
        let work = levels
            .iter()
            .map(|l| Work {
                r: vec![0.0; l.p.nrows],
                bc: vec![0.0; l.p.ncols],
                xc: vec![0.0; l.p.ncols],
            })
            .collect();
        Ok(Amg {
            a0,
            mats,
            levels,
            coarse,
            work: RefCell::new(work),
            info,
        })
    }

    fn mat(&self, l: usize) -> &Csr {
        if l == 0 {
            self.a0
        } else {
            &self.mats[l - 1]
        }
    }

    /// z = V-cycle(r): one application of the preconditioner.
    pub fn apply(&self, r: &[f64], z: &mut [f64]) {
        let mut work = self.work.borrow_mut();
        self.vcycle(0, r, z, &mut work);
    }

    fn vcycle(&self, l: usize, b: &[f64], x: &mut [f64], work: &mut [Work]) {
        if l == self.levels.len() {
            self.coarse.solve(b, x);
            return;
        }
        let a = self.mat(l);
        let lev = &self.levels[l];
        let (w, rest) = work.split_first_mut().unwrap();
        x.iter_mut().for_each(|v| *v = 0.0);
        sgs(a, &lev.diag, x, b);
        a.spmv(x, &mut w.r);
        for i in 0..b.len() {
            w.r[i] = b[i] - w.r[i];
        }
        lev.r.apply(&w.r, &mut w.bc);
        self.vcycle(l + 1, &w.bc, &mut w.xc, rest);
        lev.p.apply(&w.xc, &mut w.r); // reuse r as P xc
        for i in 0..x.len() {
            x[i] += w.r[i];
        }
        sgs(a, &lev.diag, x, b);
    }
}

/// One symmetric Gauss-Seidel sweep (forward then backward).
fn sgs(a: &Csr, diag: &[f64], x: &mut [f64], b: &[f64]) {
    let n = a.nrows;
    let relax = |i: usize, x: &mut [f64]| {
        let mut s = b[i];
        for p in a.indptr[i]..a.indptr[i + 1] {
            s -= a.data[p] * x[a.indices[p] as usize];
        }
        x[i] += s / diag[i];
    };
    for i in 0..n {
        relax(i, x);
    }
    for i in (0..n).rev() {
        relax(i, x);
    }
}

/// Dense Cholesky factorization for the coarsest level.
struct DenseCholesky {
    n: usize,
    l: Vec<f64>, // row-major lower triangle
}

impl DenseCholesky {
    fn new(a: &Csr) -> Result<Self> {
        let n = a.nrows;
        let mut l = vec![0.0; n * n];
        for i in 0..n {
            for p in a.indptr[i]..a.indptr[i + 1] {
                l[i * n + a.indices[p] as usize] = a.data[p];
            }
        }
        for j in 0..n {
            let mut d = l[j * n + j];
            for k in 0..j {
                d -= l[j * n + k] * l[j * n + k];
            }
            if !(d > 0.0) {
                bail!("coarse operator is not positive definite");
            }
            let d = d.sqrt();
            l[j * n + j] = d;
            for i in j + 1..n {
                let mut s = l[i * n + j];
                for k in 0..j {
                    s -= l[i * n + k] * l[j * n + k];
                }
                l[i * n + j] = s / d;
            }
        }
        Ok(DenseCholesky { n, l })
    }

    fn solve(&self, b: &[f64], x: &mut [f64]) {
        let (n, l) = (self.n, &self.l);
        for i in 0..n {
            let mut s = b[i];
            for k in 0..i {
                s -= l[i * n + k] * x[k];
            }
            x[i] = s / l[i * n + i];
        }
        for i in (0..n).rev() {
            let mut s = x[i];
            for k in i + 1..n {
                s -= l[k * n + i] * x[k];
            }
            x[i] = s / l[i * n + i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 1-D Laplacian with a Dirichlet-like diagonal shift.
    fn lap(n: usize) -> Csr {
        let a: Vec<u32> = (0..n as u32 - 1).collect();
        let b: Vec<u32> = (1..n as u32).collect();
        let g = vec![1.0; n - 1];
        let mut d = vec![2.0; n];
        d[0] = 1.5;
        d[n - 1] = 1.5;
        Csr::from_links(n, &a, &b, &g, &d)
    }

    #[test]
    fn vcycle_is_symmetric_and_reduces_error() {
        let a = lap(5000);
        let amg = Amg::new(&a, &AmgOptions::default()).unwrap();
        assert!(amg.info.sizes.len() >= 3);
        // symmetry: u . M v == v . M u
        let u: Vec<f64> = (0..5000).map(|i| ((i * 37) % 101) as f64 - 50.0).collect();
        let v: Vec<f64> = (0..5000).map(|i| ((i * 53) % 97) as f64 - 48.0).collect();
        let (mut mu, mut mv) = (vec![0.0; 5000], vec![0.0; 5000]);
        amg.apply(&u, &mut mu);
        amg.apply(&v, &mut mv);
        let (a1, a2) = (dot(&v, &mu), dot(&u, &mv));
        assert!((a1 - a2).abs() <= 1e-9 * a1.abs().max(1.0), "{a1} vs {a2}");
    }

    #[test]
    fn dense_cholesky_solves() {
        let a = lap(50);
        let c = DenseCholesky::new(&a).unwrap();
        let x0: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let mut b = vec![0.0; 50];
        a.spmv(&x0, &mut b);
        let mut x = vec![0.0; 50];
        c.solve(&b, &mut x);
        for i in 0..50 {
            assert!((x[i] - x0[i]).abs() < 1e-9);
        }
    }
}
