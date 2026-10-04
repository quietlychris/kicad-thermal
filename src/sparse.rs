//! Minimal CSR sparse matrix: what the CG solver and the AMG hierarchy need.
//!
//! Parallel kernels split work into fixed row chunks and combine partial results in chunk
//! order, so results are bit-for-bit identical for any number of threads.

use rayon::prelude::*;

/// Rows per parallel work item.
pub const CHUNK: usize = 8192;

#[derive(Clone, Debug)]
pub struct Csr {
    pub nrows: usize,
    pub ncols: usize,
    pub indptr: Vec<usize>,
    pub indices: Vec<u32>,
    pub data: Vec<f64>,
}

impl Csr {
    pub fn nnz(&self) -> usize {
        self.data.len()
    }

    /// The SPD matrix `diag(d) - offdiag`, with off-diagonal couplings `g` between `a[e]`
    /// and `b[e]` (each link given once; it is placed in both rows). Columns sorted per row.
    pub fn from_links(n: usize, a: &[u32], b: &[u32], g: &[f64], d: &[f64]) -> Csr {
        let mut deg = vec![1usize; n]; // diagonal
        for e in 0..a.len() {
            deg[a[e] as usize] += 1;
            deg[b[e] as usize] += 1;
        }
        let mut indptr = vec![0usize; n + 1];
        for i in 0..n {
            indptr[i + 1] = indptr[i] + deg[i];
        }
        let nnz = indptr[n];
        let mut indices = vec![0u32; nnz];
        let mut data = vec![0.0f64; nnz];
        let mut fill: Vec<usize> = indptr[..n].to_vec();
        for i in 0..n {
            indices[fill[i]] = i as u32;
            data[fill[i]] = d[i];
            fill[i] += 1;
        }
        for e in 0..a.len() {
            let (ra, rb) = (a[e] as usize, b[e] as usize);
            indices[fill[ra]] = b[e];
            data[fill[ra]] = -g[e];
            fill[ra] += 1;
            indices[fill[rb]] = a[e];
            data[fill[rb]] = -g[e];
            fill[rb] += 1;
        }
        let mut m = Csr {
            nrows: n,
            ncols: n,
            indptr,
            indices,
            data,
        };
        m.sort_rows();
        m
    }

    /// Sort each row's entries by column (parallel, deterministic).
    pub fn sort_rows(&mut self) {
        let indptr = &self.indptr;
        let rows: Vec<(usize, usize)> = (0..self.nrows)
            .map(|i| (indptr[i], indptr[i + 1]))
            .collect();
        let mut pairs: Vec<(u32, f64)> = self
            .indices
            .iter()
            .copied()
            .zip(self.data.iter().copied())
            .collect();
        // split `pairs` into per-row slices and sort them in parallel
        let mut slices: Vec<&mut [(u32, f64)]> = Vec::with_capacity(self.nrows);
        let mut rest: &mut [(u32, f64)] = &mut pairs;
        for &(s, e) in &rows {
            let (head, tail) = rest.split_at_mut(e - s);
            slices.push(head);
            rest = tail;
        }
        slices
            .par_iter_mut()
            .for_each(|s| s.sort_unstable_by_key(|p| p.0));
        for (k, (c, v)) in pairs.into_iter().enumerate() {
            self.indices[k] = c;
            self.data[k] = v;
        }
    }

    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.nrows)
            .map(|i| {
                let (s, e) = (self.indptr[i], self.indptr[i + 1]);
                (s..e)
                    .find(|&p| self.indices[p] as usize == i)
                    .map_or(0.0, |p| self.data[p])
            })
            .collect()
    }

    /// y = A x
    pub fn spmv(&self, x: &[f64], y: &mut [f64]) {
        y.par_chunks_mut(CHUNK).enumerate().for_each(|(c, ys)| {
            let r0 = c * CHUNK;
            for (o, yi) in ys.iter_mut().enumerate() {
                let i = r0 + o;
                let mut s = 0.0;
                for p in self.indptr[i]..self.indptr[i + 1] {
                    s += self.data[p] * x[self.indices[p] as usize];
                }
                *yi = s;
            }
        });
    }

    /// y = A x, returning x . y (fused for CG)
    pub fn spmv_dot(&self, x: &[f64], y: &mut [f64]) -> f64 {
        let parts: Vec<f64> = y
            .par_chunks_mut(CHUNK)
            .enumerate()
            .map(|(c, ys)| {
                let r0 = c * CHUNK;
                let mut acc = 0.0;
                for (o, yi) in ys.iter_mut().enumerate() {
                    let i = r0 + o;
                    let mut s = 0.0;
                    for p in self.indptr[i]..self.indptr[i + 1] {
                        s += self.data[p] * x[self.indices[p] as usize];
                    }
                    *yi = s;
                    acc += x[i] * s;
                }
                acc
            })
            .collect();
        parts.iter().sum()
    }

    /// y = A x for a rectangular operator (P or R); same kernel, explicit name for clarity.
    pub fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.spmv(x, y)
    }

    pub fn transpose(&self) -> Csr {
        let mut cnt = vec![0usize; self.ncols + 1];
        for &c in &self.indices {
            cnt[c as usize + 1] += 1;
        }
        for i in 0..self.ncols {
            cnt[i + 1] += cnt[i];
        }
        let indptr = cnt.clone();
        let mut fill = cnt;
        let mut indices = vec![0u32; self.nnz()];
        let mut data = vec![0.0; self.nnz()];
        for i in 0..self.nrows {
            for p in self.indptr[i]..self.indptr[i + 1] {
                let c = self.indices[p] as usize;
                indices[fill[c]] = i as u32;
                data[fill[c]] = self.data[p];
                fill[c] += 1;
            }
        }
        Csr {
            nrows: self.ncols,
            ncols: self.nrows,
            indptr,
            indices,
            data,
        }
    }

    /// C = A B (Gustavson, parallel over row chunks, columns sorted per row).
    pub fn matmul(&self, b: &Csr) -> Csr {
        assert_eq!(self.ncols, b.nrows);
        let nchunks = self.nrows.div_ceil(CHUNK).max(1);
        let parts: Vec<(Vec<usize>, Vec<u32>, Vec<f64>)> = (0..nchunks)
            .into_par_iter()
            .map_init(
                || {
                    (
                        vec![0.0f64; b.ncols],
                        vec![u32::MAX; b.ncols],
                        Vec::<u32>::new(),
                    )
                },
                |(acc, mark, cols), c| {
                    let r0 = c * CHUNK;
                    let r1 = (r0 + CHUNK).min(self.nrows);
                    let mut lens = Vec::with_capacity(r1 - r0);
                    let mut idx = Vec::new();
                    let mut val = Vec::new();
                    for i in r0..r1 {
                        cols.clear();
                        for p in self.indptr[i]..self.indptr[i + 1] {
                            let k = self.indices[p] as usize;
                            let av = self.data[p];
                            for q in b.indptr[k]..b.indptr[k + 1] {
                                let j = b.indices[q];
                                if mark[j as usize] != i as u32 {
                                    mark[j as usize] = i as u32;
                                    acc[j as usize] = 0.0;
                                    cols.push(j);
                                }
                                acc[j as usize] += av * b.data[q];
                            }
                        }
                        cols.sort_unstable();
                        let mut len = 0;
                        for &j in cols.iter() {
                            let v = acc[j as usize];
                            if v != 0.0 {
                                idx.push(j);
                                val.push(v);
                                len += 1;
                            }
                        }
                        lens.push(len);
                    }
                    (lens, idx, val)
                },
            )
            .collect();
        let mut indptr = Vec::with_capacity(self.nrows + 1);
        indptr.push(0usize);
        let total: usize = parts.iter().map(|p| p.1.len()).sum();
        let mut indices = Vec::with_capacity(total);
        let mut data = Vec::with_capacity(total);
        for (lens, idx, val) in parts {
            for l in lens {
                indptr.push(indptr[indptr.len() - 1] + l);
            }
            indices.extend(idx);
            data.extend(val);
        }
        Csr {
            nrows: self.nrows,
            ncols: b.ncols,
            indptr,
            indices,
            data,
        }
    }
}

// ------------------------------------------------------------------ dense vector kernels

/// Deterministic parallel sum of f(i) over chunked index ranges.
pub fn chunked_sum(n: usize, f: impl Fn(std::ops::Range<usize>) -> f64 + Sync) -> f64 {
    let nchunks = n.div_ceil(CHUNK);
    let parts: Vec<f64> = (0..nchunks)
        .into_par_iter()
        .map(|c| f(c * CHUNK..((c + 1) * CHUNK).min(n)))
        .collect();
    parts.iter().sum()
}

pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    chunked_sum(a.len(), |r| r.map(|i| a[i] * b[i]).sum())
}

pub fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}
