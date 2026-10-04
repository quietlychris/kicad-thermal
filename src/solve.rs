//! DC conduction on one net and steady-state heat conduction on the whole board
//! (`solve_dc`, `solve_thermal` and `_links` in kicad_fdtd's thermal.py).

use crate::cg::{self, Precond, SolveStats};
use crate::grid::{CellKind, Grid, RHO_CU};
use crate::sparse::Csr;
use anyhow::{bail, Result};

/// Conductances between neighbouring cells: link endpoints (flat indices) and G (W/K or S).
/// Order matches kicad_fdtd: x links, then y, then z; C order within each axis.
pub struct Links {
    pub a: Vec<u32>,
    pub b: Vec<u32>,
    pub g: Vec<f64>,
}

pub fn links(gr: &Grid, kx: &[f64], ky: &[f64], kz: &[f64], active: Option<&[bool]>) -> Links {
    let (nx, ny, nz) = (gr.nx, gr.ny, gr.nz);
    let dx: Vec<f64> = gr.dx.iter().map(|v| v * 1e-3).collect();
    let dy: Vec<f64> = gr.dy.iter().map(|v| v * 1e-3).collect();
    let dz: Vec<f64> = gr.dz.iter().map(|v| v * 1e-3).collect();
    let cap = 3 * nx * ny * nz;
    let (mut la, mut lb, mut lg) = (
        Vec::with_capacity(cap),
        Vec::with_capacity(cap),
        Vec::with_capacity(cap),
    );
    let flat = |i: usize, j: usize, k: usize| (i * ny + j) * nz + k;
    for axis in 0..3 {
        let kk = [kx, ky, kz][axis];
        let (ni, nj, nk) = match axis {
            0 => (nx - 1, ny, nz),
            1 => (nx, ny - 1, nz),
            _ => (nx, ny, nz.saturating_sub(1)),
        };
        for i in 0..ni {
            for j in 0..nj {
                for k in 0..nk {
                    let a = flat(i, j, k);
                    let (b, area, da, db) = match axis {
                        0 => (flat(i + 1, j, k), dy[j] * dz[k], dx[i], dx[i + 1]),
                        1 => (flat(i, j + 1, k), dx[i] * dz[k], dy[j], dy[j + 1]),
                        _ => (flat(i, j, k + 1), dx[i] * dy[j], dz[k], dz[k + 1]),
                    };
                    let mut g = area / (da / (2.0 * kk[a]) + db / (2.0 * kk[b]));
                    if !g.is_finite() {
                        g = 0.0;
                    }
                    if !(g > 0.0) {
                        continue;
                    }
                    if let Some(act) = active {
                        if !(act[a] && act[b]) {
                            continue;
                        }
                    }
                    la.push(a as u32);
                    lb.push(b as u32);
                    lg.push(g);
                }
            }
        }
    }
    Links {
        a: la,
        b: lb,
        g: lg,
    }
}

/// Diagonal of the conductance matrix: sum of G over each node's links (a side, then b side).
fn link_diagonal(n: usize, a: &[u32], b: &[u32], g: &[f64]) -> Vec<f64> {
    let mut d = vec![0.0; n];
    for e in 0..a.len() {
        d[a[e] as usize] += g[e];
    }
    for e in 0..b.len() {
        d[b[e] as usize] += g[e];
    }
    d
}

pub struct DcResult {
    pub net: String,
    pub current: f64,
    pub drop_v: f64,
    pub r_ohm: f64,
    /// volts per cell, NaN off the net's connected copper
    pub v: Vec<f64>,
    /// A/mm2 per cell
    pub j: Vec<f64>,
    /// W per cell
    pub joule: Vec<f64>,
    pub joule_w: f64,
    pub jmax: f64,
    pub cells: usize,
    pub stats: SolveStats,
}

/// DC conduction on one net. `src`/`sink` = (footprint ref, pad number).
pub fn solve_dc(
    gr: &Grid,
    net: &str,
    src: (&str, &str),
    sink: (&str, &str),
    current: f64,
    precond: Precond,
) -> Result<DcResult> {
    let (nx, ny, nz) = (gr.nx, gr.ny, gr.nz);
    let n = gr.n();
    let sig = 1.0 / RHO_CU;
    // in-plane: sheet conductance sigma * t * f ; vertical: via barrels
    let mut sx = vec![0.0; n];
    let mut sz = vec![0.0; n];
    let mut active = vec![false; n];
    for (k, kd) in gr.kind.iter().enumerate() {
        if let CellKind::Copper(layer) = *kd {
            let f = gr.copper_fraction(layer, Some(net));
            for i in 0..nx {
                for jj in 0..ny {
                    let fv = f[i * ny + jj];
                    let c = (i * ny + jj) * nz + k;
                    sx[c] = sig * fv;
                    active[c] = fv > 1e-3;
                }
            }
        }
    }
    let a_via = gr.via_copper_area(Some(net));
    for i in 0..nx {
        for jj in 0..ny {
            let cell_a = gr.dx[i] * gr.dy[jj];
            for k in 0..nz {
                let c = (i * ny + jj) * nz + k;
                if a_via[c] > 0.0 {
                    sz[c] = sig * (a_via[c] / cell_a);
                    active[c] = true;
                }
            }
        }
    }
    for c in 0..n {
        if !active[c] {
            sx[c] = 0.0;
            sz[c] = 0.0;
        }
    }
    let lk = links(gr, &sx, &sx, &sz, Some(&active));

    let flat = |c: (usize, usize, usize)| (c.0 * ny + c.1) * nz + c.2;
    let src_cells: Vec<_> = gr
        .pad_cells(src.0, Some(src.1))?
        .into_iter()
        .filter(|&c| active[flat(c)])
        .collect();
    let sink_cells: Vec<_> = gr
        .pad_cells(sink.0, Some(sink.1))?
        .into_iter()
        .filter(|&c| active[flat(c)])
        .collect();
    if src_cells.is_empty() || sink_cells.is_empty() {
        bail!(
            "source or sink pad has no copper of net {} under it at this cell size",
            net
        );
    }

    // connected component of the source (BFS over the link graph)
    let mut deg = vec![0usize; n + 1];
    for e in 0..lk.a.len() {
        deg[lk.a[e] as usize + 1] += 1;
        deg[lk.b[e] as usize + 1] += 1;
    }
    for i in 0..n {
        deg[i + 1] += deg[i];
    }
    let mut adj = vec![0u32; deg[n]];
    let mut fill = deg[..n].to_vec();
    for e in 0..lk.a.len() {
        let (a, b) = (lk.a[e] as usize, lk.b[e] as usize);
        adj[fill[a]] = b as u32;
        fill[a] += 1;
        adj[fill[b]] = a as u32;
        fill[b] += 1;
    }
    let mut in_comp = vec![false; n];
    let start = flat(src_cells[0]);
    in_comp[start] = true;
    let mut stack = vec![start];
    while let Some(u) = stack.pop() {
        for &w in &adj[deg[u]..deg[u + 1]] {
            if !in_comp[w as usize] {
                in_comp[w as usize] = true;
                stack.push(w as usize);
            }
        }
    }
    drop(adj);
    if sink_cells.iter().any(|&c| !in_comp[flat(c)]) {
        bail!("source and sink are not connected through copper of {} on this grid (try a smaller --cell)", net);
    }
    let members: Vec<usize> = (0..n).filter(|&c| in_comp[c]).collect();
    let mut loc = vec![u32::MAX; n];
    for (m, &c) in members.iter().enumerate() {
        loc[c] = m as u32;
    }
    let m = members.len();
    let keep: Vec<usize> = (0..lk.a.len())
        .filter(|&e| in_comp[lk.a[e] as usize])
        .collect();
    let a: Vec<u32> = keep.iter().map(|&e| lk.a[e]).collect();
    let b: Vec<u32> = keep.iter().map(|&e| lk.b[e]).collect();
    let g: Vec<f64> = keep.iter().map(|&e| lk.g[e]).collect();
    let ra: Vec<u32> = a.iter().map(|&c| loc[c as usize]).collect();
    let rb: Vec<u32> = b.iter().map(|&c| loc[c as usize]).collect();
    let mut diag = link_diagonal(m, &ra, &rb, &g);
    let mut rhs = vec![0.0; m];
    for &c in &src_cells {
        rhs[loc[flat(c)] as usize] += current / src_cells.len() as f64;
    }
    for &c in &sink_cells {
        rhs[loc[flat(c)] as usize] -= current / sink_cells.len() as f64;
    }
    let gnd = loc[flat(sink_cells[0])] as usize;
    let dmax = diag.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    diag[gnd] += 1e3 * dmax; // pin ~0 V
    let mat = Csr::from_links(m, &ra, &rb, &g, &diag);
    let (x, stats) = cg::solve(&mat, &rhs, cg::DEFAULT_TOL, cg::DEFAULT_MAXITER, precond)?;

    let mut v = vec![f64::NAN; n];
    for (mi, &c) in members.iter().enumerate() {
        v[c] = x[mi] - x[gnd];
    }
    let mut joule = vec![0.0; n];
    let mut jd = vec![0.0f64; n];
    let mut heat = Vec::with_capacity(a.len());
    for e in 0..a.len() {
        let d = x[ra[e] as usize] - x[rb[e] as usize];
        heat.push(0.5 * g[e] * (d * d));
    }
    for e in 0..a.len() {
        joule[a[e] as usize] += heat[e];
    }
    for e in 0..b.len() {
        joule[b[e] as usize] += heat[e];
    }
    // current density: in-plane links only, J = I / (t * f * width), width = transverse cell size
    for e in 0..a.len() {
        let (ca, cb) = (a[e] as usize, b[e] as usize);
        let (ka, kb) = (ca % nz, cb % nz);
        let (ia, ja) = ((ca / nz) / ny, (ca / nz) % ny);
        let ib = (cb / nz) / ny;
        let i_link = g[e] * (x[ra[e] as usize] - x[rb[e] as usize]);
        let t = gr.dz[ka] * 1e-3;
        let f = sx[ca].max(sx[cb]) / sig;
        let width = if ia != ib {
            gr.dy[ja.min(ny - 1)]
        } else {
            gr.dx[ia.min(nx - 1)]
        } * 1e-3;
        let jl = if ka == kb && f > 0.0 {
            i_link.abs() / (t * f * width)
        } else {
            0.0
        } / 1e6;
        if jl > jd[ca] {
            jd[ca] = jl;
        }
        if jl > jd[cb] {
            jd[cb] = jl;
        }
    }
    let vs = src_cells.iter().map(|&c| v[flat(c)]).sum::<f64>() / src_cells.len() as f64;
    Ok(DcResult {
        net: net.to_string(),
        current,
        drop_v: vs,
        r_ohm: if current != 0.0 { vs / current } else { 0.0 },
        joule_w: joule.iter().sum(),
        jmax: jd.iter().cloned().fold(0.0, f64::max),
        v,
        j: jd,
        joule,
        cells: m,
        stats,
    })
}

pub struct ThermalResult {
    /// temperature per cell, C
    pub t: Vec<f64>,
    /// heat input per cell, W
    pub q: Vec<f64>,
    pub h: f64,
    pub tamb: f64,
    pub p_total: f64,
    pub tmax: f64,
    pub stats: SolveStats,
}

/// Steady-state temperature. `powers`: (footprint ref, W) in order; `extra_q`: per-cell W
/// (e.g. DC Joule heat).
pub fn solve_thermal(
    gr: &Grid,
    powers: &[(String, f64)],
    h: f64,
    tamb: f64,
    extra_q: Option<&[f64]>,
    precond: Precond,
) -> Result<ThermalResult> {
    let (nx, ny, nz) = (gr.nx, gr.ny, gr.nz);
    let n = gr.n();
    let lk = links(gr, &gr.kx, &gr.kx, &gr.kz, None);
    let mut diag = link_diagonal(n, &lk.a, &lk.b, &lk.g);
    let mut q = vec![0.0; n];
    for (reference, p) in powers {
        let cells = gr.pad_cells(reference, None)?;
        if cells.is_empty() {
            bail!("no pads for footprint {}", reference);
        }
        for &(i, j, k) in &cells {
            q[(i * ny + j) * nz + k] += p / cells.len() as f64;
        }
    }
    if let Some(e) = extra_q {
        for c in 0..n {
            q[c] += e[c];
        }
    }
    // convection on the top face of the top cells and the bottom face of the bottom cells
    let mut ha = vec![0.0; n];
    for i in 0..nx {
        for j in 0..ny {
            let face_a = (gr.dx[i] * gr.dy[j]) * 1e-6;
            ha[(i * ny + j) * nz] += h * face_a;
            ha[(i * ny + j) * nz + nz - 1] += h * face_a;
        }
    }
    let mut rhs = vec![0.0; n];
    for c in 0..n {
        diag[c] += ha[c];
        rhs[c] = q[c] + ha[c] * tamb;
    }
    let mat = Csr::from_links(n, &lk.a, &lk.b, &lk.g, &diag);
    drop(lk);
    let (t, stats) = cg::solve(&mat, &rhs, cg::DEFAULT_TOL, cg::DEFAULT_MAXITER, precond)?;
    let tmax = t.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    Ok(ThermalResult {
        p_total: q.iter().sum(),
        t,
        q,
        h,
        tamb,
        tmax,
        stats,
    })
}
