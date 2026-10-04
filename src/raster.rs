//! Point-in-polygon-union sampling along horizontal scanlines.
//!
//! kicad_fdtd samples copper with `shapely.contains_xy(unary_union(polys), X, Y)`: a point
//! counts when it lies strictly inside the union. Here each polygon (outline plus holes,
//! even-odd) is cut by the scanline into open interior intervals; the intervals of all
//! polygons are merged, with touching intervals joined (a shared edge is interior to the
//! union), and a sample counts when it lies strictly inside a merged interval. A polygon's
//! interior is its outline (even-odd) minus the union of its holes.
//!
//! Differences from shapely, both measure-zero for real boards: samples exactly on a
//! polygon edge may be classified differently than shapely's exact predicates would, and
//! a self-intersecting outline is filled even-odd (shapely's `buffer(0)` repair can differ).

use crate::geometry::Polygon;

/// Even-odd interior intervals of one ring on the line `y`, appended to `out`.
fn ring_intervals(ring: &[(f64, f64)], y: f64, xs: &mut Vec<f64>, out: &mut Vec<(f64, f64)>) {
    xs.clear();
    let n = ring.len();
    for i in 0..n {
        let (x1, y1) = ring[i];
        let (x2, y2) = ring[(i + 1) % n];
        if (y1 > y) != (y2 > y) {
            xs.push(x1 + (y - y1) * (x2 - x1) / (y2 - y1));
        }
    }
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    for pair in xs.chunks_exact(2) {
        if pair[1] > pair[0] {
            out.push((pair[0], pair[1]));
        }
    }
}

/// Sort and merge open intervals in place; touching intervals join (their shared end is interior).
fn merge(iv: &mut Vec<(f64, f64)>) {
    if iv.len() < 2 {
        return;
    }
    iv.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mut w = 0;
    for r in 1..iv.len() {
        let (a, b) = iv[r];
        if a <= iv[w].1 {
            if b > iv[w].1 {
                iv[w].1 = b;
            }
        } else {
            w += 1;
            iv[w] = (a, b);
        }
    }
    iv.truncate(w + 1);
}

/// Open interior intervals of the union of `polys` on the horizontal line `y`, merged and sorted.
///
/// Each polygon is its outline minus the union of its holes, which is what shapely's
/// `buffer(0)` repair gives for the common invalid case of overlapping holes.
pub fn union_intervals(polys: &[&Polygon], y: f64, xs: &mut Vec<f64>, out: &mut Vec<(f64, f64)>) {
    out.clear();
    let mut outer = Vec::new();
    let mut holes = Vec::new();
    for p in polys {
        if y < p.ymin || y > p.ymax {
            continue;
        }
        outer.clear();
        holes.clear();
        ring_intervals(&p.rings[0], y, xs, &mut outer);
        if outer.is_empty() {
            continue;
        }
        for h in &p.rings[1..] {
            ring_intervals(h, y, xs, &mut holes);
        }
        merge(&mut holes);
        // outer minus holes (both open intervals; hole boundaries are not interior)
        let mut hi = 0;
        for &(a, b) in &outer {
            let mut start = a;
            while hi < holes.len() && holes[hi].1 <= start {
                hi += 1;
            }
            let mut k = hi;
            while k < holes.len() && holes[k].0 < b {
                let (c, d) = holes[k];
                if c > start {
                    out.push((start, c));
                }
                if d > start {
                    start = d;
                }
                k += 1;
            }
            if start < b {
                out.push((start, b));
            }
        }
    }
    merge(out);
}

/// For sorted sample abscissae `xs`, call `hit(index)` for each sample strictly inside an interval.
pub fn for_each_inside(intervals: &[(f64, f64)], xs: &[f64], mut hit: impl FnMut(usize)) {
    for &(a, b) in intervals {
        let lo = xs.partition_point(|&x| x <= a);
        let hi = xs.partition_point(|&x| x < b);
        for s in lo..hi {
            hit(s);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rect(xa: f64, ya: f64, xb: f64, yb: f64) -> Polygon {
        Polygon::new(vec![vec![(xa, ya), (xb, ya), (xb, yb), (xa, yb)]])
    }

    #[test]
    fn touching_rectangles_merge_and_holes_cut() {
        let a = rect(0.0, 0.0, 1.0, 1.0);
        let b = rect(1.0, 0.0, 2.0, 1.0);
        let mut c = rect(3.0, 0.0, 6.0, 1.0);
        c.rings
            .push(vec![(4.0, 0.2), (5.0, 0.2), (5.0, 0.8), (4.0, 0.8)]);
        let (mut xs, mut out) = (Vec::new(), Vec::new());
        union_intervals(&[&a, &b, &c], 0.5, &mut xs, &mut out);
        assert_eq!(out, vec![(0.0, 2.0), (3.0, 4.0), (5.0, 6.0)]);
        let samples = [0.0, 1.0, 2.0, 3.5, 4.5, 5.5];
        let mut hits = Vec::new();
        for_each_inside(&out, &samples, |s| hits.push(s));
        assert_eq!(hits, vec![1, 3, 5]); // 0.0 and 2.0 are on the boundary; 4.5 is in the hole
    }

    #[test]
    fn overlapping_holes_are_unioned_not_even_odd() {
        let mut p = rect(0.0, 0.0, 10.0, 1.0);
        p.rings
            .push(vec![(2.0, 0.0), (5.0, 0.0), (5.0, 1.0), (2.0, 1.0)]);
        p.rings
            .push(vec![(4.0, 0.0), (7.0, 0.0), (7.0, 1.0), (4.0, 1.0)]);
        let (mut xs, mut out) = (Vec::new(), Vec::new());
        union_intervals(&[&p], 0.5, &mut xs, &mut out);
        assert_eq!(out, vec![(0.0, 2.0), (7.0, 10.0)]);
    }
}
