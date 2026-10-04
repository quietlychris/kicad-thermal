//! Plain-text report and result JSON, in the same format as kicad_fdtd's thermal CLI.

use crate::grid::Grid;
use crate::solve::{DcResult, ThermalResult};
use anyhow::Result;
use serde_json::{json, Map, Value};

pub fn report(gr: &Grid, th: &ThermalResult, dc: Option<&DcResult>) -> Result<String> {
    let (ny, nz) = (gr.ny, gr.nz);
    let t = &th.t;
    let mut lines = vec![format!(
        "thermal: {:.2} W in, h {:.0} W/m2K both faces, T_amb {:.0} C, grid {} x {} x {} cells of {:.2} mm",
        th.p_total, th.h, th.tamb, gr.nx, gr.ny, gr.nz, gr.cell
    )];
    // first maximum in C order, like numpy.argmax
    let mut imax = 0;
    for c in 1..t.len() {
        if t[c] > t[imax] {
            imax = c;
        }
    }
    let (i, j, k) = ((imax / nz) / ny, (imax / nz) % ny, imax % nz);
    lines.push(format!(
        "max {:.1} C (+{:.1} K) at ({:.1}, {:.1}) mm on {}",
        t[imax],
        t[imax] - th.tamb,
        gr.cx[i],
        gr.cy[j],
        gr.layer_name_of_k(k)
    ));
    for (l, layer) in gr.bd.layers.iter().enumerate() {
        let k = gr.k_cu[l];
        let (mut mx, mut sum, mut cnt) = (f64::NEG_INFINITY, 0.0, 0usize);
        for i in 0..gr.nx {
            for j in 0..ny {
                let v = t[(i * ny + j) * nz + k];
                mx = mx.max(v);
                if gr.inside[i * ny + j] {
                    sum += v;
                    cnt += 1;
                }
            }
        }
        lines.push(format!(
            "   {:<7} max {:.1} C  mean {:.1} C",
            layer,
            mx,
            sum / cnt as f64
        ));
    }
    let mut rows: Vec<(f64, String, String)> = Vec::new();
    for fp in &gr.bd.footprints {
        let cells = gr.pad_cells(&fp.reference, None)?;
        if !cells.is_empty() {
            let tm = cells
                .iter()
                .map(|&(i, j, k)| t[(i * ny + j) * nz + k])
                .fold(f64::NEG_INFINITY, f64::max);
            rows.push((tm, fp.reference.clone(), fp.value.clone()));
        }
    }
    rows.sort_by(|a, b| {
        b.0.partial_cmp(&a.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| b.1.cmp(&a.1))
            .then_with(|| b.2.cmp(&a.2))
    });
    lines.push(format!(
        "hottest footprints: {}",
        rows.iter()
            .take(8)
            .map(|(t, r, v)| format!("{} ({}) {:.1} C", r, v, t))
            .collect::<Vec<_>>()
            .join(", ")
    ));
    if let Some(dc) = dc {
        lines.push(format!(
            "DC {}: {:.2} A from source to sink, drop {:.1} mV ({:.1} mohm), Joule {:.3} W, max current density {:.1} A/mm2 over {} copper cells",
            dc.net,
            dc.current,
            1e3 * dc.drop_v,
            1e3 * dc.r_ohm,
            dc.joule_w,
            dc.jmax,
            dc.cells
        ));
    }
    Ok(lines.join("\n"))
}

/// One copper layer of a per-cell field as nested [nx][ny] lists; NaN becomes `nan_as`.
fn layer_field(gr: &Grid, f: &[f64], k: usize, nan_as: f64) -> Value {
    let (ny, nz) = (gr.ny, gr.nz);
    Value::Array(
        (0..gr.nx)
            .map(|i| {
                Value::Array(
                    (0..ny)
                        .map(|j| {
                            let v = f[(i * ny + j) * nz + k];
                            json!(if v.is_nan() { nan_as } else { v })
                        })
                        .collect(),
                )
            })
            .collect(),
    )
}

fn per_layer(gr: &Grid, f: &[f64], nan_as: f64) -> Value {
    let mut m = Map::new();
    for (l, layer) in gr.bd.layers.iter().enumerate() {
        m.insert(layer.clone(), layer_field(gr, f, gr.k_cu[l], nan_as));
    }
    Value::Object(m)
}

/// The JSON kicad_fdtd writes to out/thermal_<board>.json.
pub fn result_json(gr: &Grid, th: &ThermalResult, dc: Option<&DcResult>) -> Result<Value> {
    let mut out = Map::new();
    out.insert("cell".into(), json!(gr.cell));
    out.insert("x".into(), json!(gr.x));
    out.insert("y".into(), json!(gr.y));
    out.insert("layers".into(), json!(gr.bd.layers));
    out.insert("h".into(), json!(th.h));
    out.insert("tamb".into(), json!(th.tamb));
    out.insert("P_total".into(), json!(th.p_total));
    out.insert("Tmax".into(), json!(th.tmax));
    out.insert("report".into(), json!(report(gr, th, dc)?));
    out.insert("T".into(), per_layer(gr, &th.t, f64::NAN));
    if let Some(dc) = dc {
        let mut d = Map::new();
        d.insert("net".into(), json!(dc.net));
        d.insert("current".into(), json!(dc.current));
        d.insert("drop_V".into(), json!(dc.drop_v));
        d.insert("R_ohm".into(), json!(dc.r_ohm));
        d.insert("joule_W".into(), json!(dc.joule_w));
        d.insert("Jmax".into(), json!(dc.jmax));
        d.insert("V".into(), per_layer(gr, &dc.v, -1.0));
        d.insert("J".into(), per_layer(gr, &dc.j, f64::NAN));
        out.insert("dc".into(), Value::Object(d));
    }
    Ok(Value::Object(out))
}
