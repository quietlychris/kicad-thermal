//! Whole-board rectilinear grid and material conductivities (`thermal.Grid` in kicad_fdtd).
//!
//! Cells are indexed (i, j, k) with x, y in plane and k running from the top copper
//! downwards; flat index = (i * ny + j) * nz + k, the same C order numpy uses.

use crate::geometry::Board;
use crate::raster::{for_each_inside, union_intervals};
use anyhow::{bail, Result};
use rayon::prelude::*;

pub const K_CU: f64 = 385.0; // W/mK
pub const K_FR4_XY: f64 = 0.8;
pub const K_FR4_Z: f64 = 0.3;
pub const K_AIR: f64 = 0.03;
pub const RHO_CU: f64 = 1.68e-8; // ohm m
pub const VIA_WALL: f64 = 0.025; // mm plating
pub const DZ_MAX: f64 = 0.4; // mm, max dielectric cell height

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CellKind {
    /// a copper layer, by index into `Board::layers`
    Copper(usize),
    /// a dielectric cell of dielectric `i` (between copper layers i and i+1)
    Dielectric(usize),
}

pub struct Grid<'b> {
    pub bd: &'b Board,
    pub cell: f64,
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    /// cell edges (nx+1, ny+1, nz+1); z decreasing
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    pub dx: Vec<f64>,
    pub dy: Vec<f64>,
    pub dz: Vec<f64>,
    pub cx: Vec<f64>,
    pub cy: Vec<f64>,
    pub kind: Vec<CellKind>,
    /// k index of each copper layer, indexed like `Board::layers`
    pub k_cu: Vec<usize>,
    /// in-plane (kx = ky) and through-plane conductivity per cell, W/mK
    pub kx: Vec<f64>,
    pub kz: Vec<f64>,
    /// copper fraction of every copper layer (all nets), nx*ny per layer
    pub cu_frac: Vec<Vec<f64>>,
    /// cell centre inside the board outline, nx*ny
    pub inside: Vec<bool>,
}

/// numpy.linspace(start, stop, num), bit for bit.
pub fn linspace(start: f64, stop: f64, num: usize) -> Vec<f64> {
    if num == 1 {
        return vec![start];
    }
    let step = (stop - start) / (num - 1) as f64;
    let mut v: Vec<f64> = (0..num).map(|i| i as f64 * step + start).collect();
    v[num - 1] = stop;
    v
}

/// numpy.searchsorted(a, v) (side="left"): the number of entries strictly below v.
pub fn searchsorted_left(a: &[f64], v: f64) -> usize {
    a.partition_point(|&e| e < v)
}

impl<'b> Grid<'b> {
    pub fn new(bd: &'b Board, cell: f64) -> Result<Grid<'b>> {
        if !(cell > 0.0) {
            bail!("cell size must be positive");
        }
        let nxl = 2usize.max((bd.w / cell).ceil() as usize + 1);
        let nyl = 2usize.max((bd.h / cell).ceil() as usize + 1);
        let x = linspace(0.0, bd.w, nxl);
        let y = linspace(-bd.h, 0.0, nyl);
        let mut zs = vec![0.0f64];
        let mut kind = Vec::new();
        let nl = bd.layers.len();
        for i in 0..nl {
            let t = bd.cu_thickness[i];
            zs.push(zs[zs.len() - 1] - t);
            kind.push(CellKind::Copper(i));
            if i < nl - 1 {
                let d = bd.diel[i];
                let span = d.z_top - d.z_bot;
                let n = 1usize.max((span / DZ_MAX).ceil() as usize);
                for _ in 0..n {
                    zs.push(zs[zs.len() - 1] - span / n as f64);
                    kind.push(CellKind::Dielectric(i));
                }
            }
        }
        let (nx, ny, nz) = (nxl - 1, nyl - 1, zs.len() - 1);
        let diff = |v: &[f64]| v.windows(2).map(|w| w[1] - w[0]).collect::<Vec<f64>>();
        let dx = diff(&x);
        let dy = diff(&y);
        let dz: Vec<f64> = zs.windows(2).map(|w| -(w[1] - w[0])).collect();
        let mid = |v: &[f64]| {
            v.windows(2)
                .map(|w| (w[0] + w[1]) / 2.0)
                .collect::<Vec<f64>>()
        };
        let cx = mid(&x);
        let cy = mid(&y);
        let mut k_cu = vec![0usize; nl];
        for (k, kd) in kind.iter().enumerate() {
            if let CellKind::Copper(l) = kd {
                k_cu[*l] = k;
            }
        }
        let mut g = Grid {
            bd,
            cell,
            nx,
            ny,
            nz,
            x,
            y,
            z: zs,
            dx,
            dy,
            dz,
            cx,
            cy,
            kind,
            k_cu,
            kx: Vec::new(),
            kz: Vec::new(),
            cu_frac: Vec::new(),
            inside: Vec::new(),
        };
        g.inside = g.inside_outline();
        g.build_materials();
        Ok(g)
    }

    #[inline]
    pub fn n(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    #[inline]
    pub fn flat(&self, i: usize, j: usize, k: usize) -> usize {
        (i * self.ny + j) * self.nz + k
    }

    /// Fraction of each cell covered by copper on one layer (all nets, or one net), sampled
    /// 3 x 3 per cell. Returns nx*ny values, index i*ny + j.
    pub fn copper_fraction(&self, layer: usize, net: Option<&str>) -> Vec<f64> {
        let polys: Vec<_> = self.bd.polys[layer]
            .iter()
            .filter(|p| net.map_or(true, |n| p.net == n))
            .map(|p| &p.poly)
            .collect();
        let (nx, ny) = (self.nx, self.ny);
        if polys.is_empty() {
            return vec![0.0; nx * ny];
        }
        let o = [-1.0 / 3.0, 0.0, 1.0 / 3.0];
        let xs: Vec<f64> = (0..nx)
            .flat_map(|i| o.iter().map(move |&ox| (i, ox)))
            .map(|(i, ox)| self.cx[i] + ox * self.dx[i])
            .collect();
        // one column (all i) of counts per j
        let cols: Vec<Vec<u8>> = (0..ny)
            .into_par_iter()
            .map_init(
                || (Vec::new(), Vec::new()),
                |(cross, ivs), j| {
                    let mut counts = vec![0u8; nx];
                    for &oy in &o {
                        let yv = self.cy[j] + oy * self.dy[j];
                        union_intervals(&polys, yv, cross, ivs);
                        for_each_inside(ivs, &xs, |s| counts[s / 3] += 1);
                    }
                    counts
                },
            )
            .collect();
        let mut f = vec![0.0; nx * ny];
        for (j, col) in cols.iter().enumerate() {
            for i in 0..nx {
                f[i * ny + j] = col[i] as f64 / 9.0;
            }
        }
        f
    }

    fn inside_outline(&self) -> Vec<bool> {
        let (nx, ny) = (self.nx, self.ny);
        let poly = [&self.bd.outline];
        let mut inside = vec![false; nx * ny];
        let (mut cross, mut ivs) = (Vec::new(), Vec::new());
        for j in 0..ny {
            union_intervals(&poly, self.cy[j], &mut cross, &mut ivs);
            for_each_inside(&ivs, &self.cx, |i| inside[i * ny + j] = true);
        }
        inside
    }

    /// Copper cross-section of via barrels per (i, j, k) cell, mm2 (all nets, or one net).
    pub fn via_copper_area(&self, net: Option<&str>) -> Vec<f64> {
        let mut a = vec![0.0; self.n()];
        for v in &self.bd.vias {
            if let Some(n) = net {
                if v.net != n {
                    continue;
                }
            }
            let i = searchsorted_left(&self.x, v.x) as isize - 1;
            let j = searchsorted_left(&self.y, v.y) as isize - 1;
            if !(0 <= i && (i as usize) < self.nx && 0 <= j && (j as usize) < self.ny) {
                continue;
            }
            let (i, j) = (i as usize, j as usize);
            // barrel = plated hole wall: drill radius + plating
            let ri = if v.drill != 0.0 {
                v.drill / 2.0
            } else {
                (v.r - VIA_WALL).max(0.0)
            };
            let ro = ri + VIA_WALL;
            let area = std::f64::consts::PI * (ro * ro - ri * ri);
            let k_of = |name: &Option<String>, default: usize| {
                name.as_deref()
                    .and_then(|l| self.bd.layer_index(l))
                    .map_or(default, |l| self.k_cu[l])
            };
            let kt = k_of(&v.top, 0);
            let kb = k_of(&v.bottom, self.nz - 1);
            for k in kt.min(kb)..=kt.max(kb) {
                let f = self.flat(i, j, k);
                a[f] += area;
            }
        }
        a
    }

    fn build_materials(&mut self) {
        let n = self.n();
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        let mut kx = vec![K_FR4_XY; n];
        let mut kz = vec![K_FR4_Z; n];
        let mut cu_frac = Vec::with_capacity(self.bd.layers.len());
        for (k, kd) in self.kind.clone().into_iter().enumerate() {
            if let CellKind::Copper(layer) = kd {
                let f = self.copper_fraction(layer, None);
                for i in 0..nx {
                    for j in 0..ny {
                        let fv = f[i * ny + j];
                        let c = (i * ny + j) * nz + k;
                        kx[c] = fv * K_CU + (1.0 - fv) * K_FR4_XY;
                        kz[c] = fv * K_CU + (1.0 - fv) * K_FR4_Z;
                    }
                }
                cu_frac.push(f);
            }
        }
        let a = self.via_copper_area(None);
        for i in 0..nx {
            for j in 0..ny {
                let cell_a = self.dx[i] * self.dy[j];
                let out = !self.inside[i * ny + j];
                for k in 0..nz {
                    let c = (i * ny + j) * nz + k;
                    kz[c] += a[c] / cell_a * K_CU;
                    if out {
                        kx[c] = K_AIR;
                        kz[c] = K_AIR;
                    }
                }
            }
        }
        self.kx = kx;
        self.kz = kz;
        self.cu_frac = cu_frac;
    }

    pub fn cell_of(&self, x: f64, y: f64) -> (usize, usize) {
        let clip = |v: isize, hi: usize| v.clamp(0, hi as isize - 1) as usize;
        (
            clip(searchsorted_left(&self.x, x) as isize - 1, self.nx),
            clip(searchsorted_left(&self.y, y) as isize - 1, self.ny),
        )
    }

    /// (i, j, k) of the cell under each pad of a footprint (or one pad), on its outer copper layer.
    pub fn pad_cells(
        &self,
        reference: &str,
        num: Option<&str>,
    ) -> Result<Vec<(usize, usize, usize)>> {
        let bd = self.bd;
        let first = &bd.layers[0];
        let last = &bd.layers[bd.layers.len() - 1];
        let mut out = Vec::new();
        for p in &bd.pads {
            if p.reference != reference || num.map_or(false, |n| p.num != n) || p.layers.is_empty()
            {
                continue; // pads with no copper layers: NPTH holes, paste-only
            }
            let lay = if p.layers.contains(first) {
                first
            } else if p.layers.contains(last) {
                last
            } else {
                &p.layers[0]
            };
            let Some(li) = bd.layer_index(lay) else {
                bail!(
                    "pad {}.{} is on unknown copper layer {}",
                    p.reference,
                    p.num,
                    lay
                );
            };
            let (i, j) = self.cell_of(p.x, p.y);
            out.push((i, j, self.k_cu[li]));
        }
        Ok(out)
    }

    /// Name of the layer a k index belongs to, or "dielectric".
    pub fn layer_name_of_k(&self, k: usize) -> &str {
        match self.kind[k] {
            CellKind::Copper(l) => &self.bd.layers[l],
            CellKind::Dielectric(_) => "dielectric",
        }
    }
}
