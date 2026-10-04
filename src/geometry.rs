//! Board geometry from kicad_fdtd's exporter JSON (`kicad-fdtd export`).
//!
//! Mirrors `kicad_fdtd/geometry.py` for everything the thermal/DC solver uses.
//! Model frame: x right, y UP (KiCad's y is flipped), z = 0 at the top copper,
//! origin at the board's bounding-box corner. All lengths in mm.

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Deserializer};
use std::collections::HashMap;
use std::path::Path;

pub const DEFAULT_ER: f64 = 4.5;
pub const DEFAULT_TAND: f64 = 0.02;
pub const DEFAULT_TOTAL: f64 = 1.6;
pub const DEFAULT_CU: f64 = 0.035;

// ------------------------------------------------------------------ raw JSON

#[derive(Deserialize)]
struct RawGeometry {
    #[serde(default)]
    board: Option<String>,
    #[serde(default)]
    kicad: Option<String>,
    layers: Vec<String>,
    bbox: [f64; 4],
    #[serde(default)]
    outline: Option<Vec<[f64; 2]>>,
    #[serde(default)]
    stackup: Option<Vec<RawStack>>,
    #[serde(default)]
    copper: HashMap<String, Vec<RawPoly>>,
    #[serde(default)]
    vias: Vec<RawVia>,
    #[serde(default)]
    pads: Vec<RawPad>,
    #[serde(default)]
    footprints: Vec<RawFootprint>,
}

#[derive(Deserialize)]
struct RawStack {
    #[serde(rename = "type")]
    typ: String,
    #[serde(default)]
    thickness_mm: Option<f64>,
    #[serde(default)]
    er: Option<f64>,
    #[serde(default)]
    tand: Option<f64>,
}

#[derive(Deserialize)]
struct RawPoly {
    #[serde(default, deserialize_with = "null_to_empty")]
    net: String,
    #[serde(default, deserialize_with = "null_to_empty")]
    kind: String,
    outline: Vec<[f64; 2]>,
    #[serde(default)]
    holes: Vec<Vec<[f64; 2]>>,
}

#[derive(Deserialize)]
struct RawVia {
    x: f64,
    y: f64,
    #[serde(default)]
    drill: Option<f64>,
    dia: f64,
    #[serde(default)]
    top: Option<String>,
    #[serde(default)]
    bottom: Option<String>,
    #[serde(default, deserialize_with = "null_to_empty")]
    net: String,
    #[serde(default, deserialize_with = "null_to_empty")]
    kind: String,
}

#[derive(Deserialize)]
struct RawPad {
    #[serde(rename = "ref")]
    reference: String,
    #[serde(deserialize_with = "string_or_number")]
    num: String,
    #[serde(default, deserialize_with = "null_to_empty")]
    net: String,
    x: f64,
    y: f64,
    #[serde(default)]
    layers: Vec<String>,
}

#[derive(Deserialize)]
struct RawFootprint {
    #[serde(rename = "ref")]
    reference: String,
    #[serde(default, deserialize_with = "null_to_empty")]
    value: String,
}

fn null_to_empty<'de, D: Deserializer<'de>>(d: D) -> std::result::Result<String, D::Error> {
    Ok(Option::<String>::deserialize(d)?.unwrap_or_default())
}

/// The exporter writes pad numbers as strings; accept numbers too (hand-written files).
fn string_or_number<'de, D: Deserializer<'de>>(d: D) -> std::result::Result<String, D::Error> {
    Ok(match serde_json::Value::deserialize(d)? {
        serde_json::Value::String(s) => s,
        serde_json::Value::Null => String::new(),
        v => v.to_string(),
    })
}

// ------------------------------------------------------------------ model types

/// A polygon with holes, in the model frame. Ring 0 is the outline, the rest are holes.
#[derive(Clone, Debug)]
pub struct Polygon {
    pub rings: Vec<Vec<(f64, f64)>>,
    pub ymin: f64,
    pub ymax: f64,
}

impl Polygon {
    pub fn new(rings: Vec<Vec<(f64, f64)>>) -> Self {
        let (mut ymin, mut ymax) = (f64::INFINITY, f64::NEG_INFINITY);
        for &(_, y) in rings.iter().flatten() {
            ymin = ymin.min(y);
            ymax = ymax.max(y);
        }
        Polygon { rings, ymin, ymax }
    }

    /// Shoelace area of the outline minus the holes (used only to drop empty polygons).
    pub fn area(&self) -> f64 {
        let ring_area = |r: &Vec<(f64, f64)>| {
            let n = r.len();
            (0..n)
                .map(|i| {
                    let (x1, y1) = r[i];
                    let (x2, y2) = r[(i + 1) % n];
                    x1 * y2 - x2 * y1
                })
                .sum::<f64>()
                .abs()
                / 2.0
        };
        let mut a = 0.0;
        for (i, r) in self.rings.iter().enumerate() {
            a += if i == 0 { ring_area(r) } else { -ring_area(r) };
        }
        a
    }
}

#[derive(Clone, Debug)]
pub struct CopperPoly {
    pub net: String,
    pub kind: String,
    pub poly: Polygon,
}

#[derive(Clone, Debug)]
pub struct Pad {
    pub reference: String,
    pub num: String,
    pub net: String,
    /// model-frame position
    pub x: f64,
    pub y: f64,
    pub layers: Vec<String>,
}

#[derive(Clone, Debug)]
pub struct Via {
    /// model-frame position
    pub x: f64,
    pub y: f64,
    /// pad (or ring) radius = dia / 2
    pub r: f64,
    pub drill: f64,
    pub top: Option<String>,
    pub bottom: Option<String>,
    pub net: String,
    pub kind: String,
}

#[derive(Clone, Debug)]
pub struct Footprint {
    pub reference: String,
    pub value: String,
}

#[derive(Clone, Copy, Debug)]
pub struct Dielectric {
    pub z_top: f64,
    pub z_bot: f64,
    pub er: f64,
    pub tand: f64,
}

/// A board loaded from the exporter's geometry JSON.
pub struct Board {
    pub board: Option<String>,
    pub kicad: Option<String>,
    pub layers: Vec<String>,
    pub x0: f64,
    pub y0: f64,
    pub w: f64,
    pub h: f64,
    /// copper thickness per layer (mm), indexed like `layers`
    pub cu_thickness: Vec<f64>,
    /// z of each copper layer (mm), indexed like `layers`
    pub z: Vec<f64>,
    pub diel: Vec<Dielectric>,
    pub z_bot: f64,
    pub stackup_source: String,
    pub outline: Polygon,
    /// copper polygons per layer, indexed like `layers`, in file order
    pub polys: Vec<Vec<CopperPoly>>,
    pub pads: Vec<Pad>,
    pub vias: Vec<Via>,
    pub footprints: Vec<Footprint>,
}

impl Board {
    pub fn load(path: impl AsRef<Path>) -> Result<Board> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("reading geometry {}", path.display()))?;
        Board::from_json_str(&text).with_context(|| format!("parsing geometry {}", path.display()))
    }

    pub fn from_json_str(text: &str) -> Result<Board> {
        let g: RawGeometry = serde_json::from_str(text)?;
        if g.layers.is_empty() {
            bail!("geometry has no copper layers");
        }
        let [x0, y0, w, h] = g.bbox;
        let n = g.layers.len();

        // --- stackup (geometry.py: Board._build_stackup)
        let stack = g.stackup.unwrap_or_default();
        let cu: Vec<&RawStack> = stack.iter().filter(|s| s.typ == "copper").collect();
        let dl: Vec<&RawStack> = stack
            .iter()
            .filter(|s| s.typ == "prepreg" || s.typ == "core")
            .collect();
        let truthy = |v: Option<f64>| v.filter(|t| *t != 0.0);
        let usable = cu.len() == n
            && dl.len() == n - 1
            && dl.iter().all(|d| truthy(d.thickness_mm).is_some());
        let (stackup_source, cu_t, dl_t): (String, Vec<f64>, Vec<(f64, f64, f64)>) = if usable {
            (
                "file".to_string(),
                cu.iter()
                    .map(|c| truthy(c.thickness_mm).unwrap_or(DEFAULT_CU))
                    .collect(),
                dl.iter()
                    .map(|d| {
                        (
                            d.thickness_mm.unwrap(),
                            truthy(d.er).unwrap_or(DEFAULT_ER),
                            truthy(d.tand).unwrap_or(DEFAULT_TAND),
                        )
                    })
                    .collect(),
            )
        } else {
            let t = (DEFAULT_TOTAL - n as f64 * DEFAULT_CU) / (n.saturating_sub(1).max(1)) as f64;
            (
                format!(
                    "default (no stackup in the file: {:.1} mm FR4, er {:.1})",
                    DEFAULT_TOTAL, DEFAULT_ER
                ),
                vec![DEFAULT_CU; n],
                vec![(t, DEFAULT_ER, DEFAULT_TAND); n - 1],
            )
        };
        let mut z = Vec::with_capacity(n);
        let mut diel = Vec::with_capacity(n.saturating_sub(1));
        let mut zc = 0.0;
        for i in 0..n {
            z.push(zc);
            if i < n - 1 {
                let (t, er, tand) = dl_t[i];
                diel.push(Dielectric {
                    z_top: zc,
                    z_bot: zc - t,
                    er,
                    tand,
                });
                zc -= t;
            }
        }

        // --- frames
        let xy = |pts: &[[f64; 2]]| -> Vec<(f64, f64)> {
            pts.iter().map(|p| (p[0] - x0, -(p[1] - y0))).collect()
        };

        let outline_pts = g.outline.unwrap_or_default();
        let outline = if outline_pts.len() >= 3 {
            Polygon::new(vec![xy(&outline_pts)])
        } else {
            Polygon::new(vec![vec![(0.0, -h), (w, -h), (w, 0.0), (0.0, 0.0)]])
        };

        let mut polys = Vec::with_capacity(n);
        for layer in &g.layers {
            let mut out = Vec::new();
            if let Some(list) = g.copper.get(layer) {
                for p in list {
                    if p.outline.len() < 3 {
                        continue;
                    }
                    let mut rings = vec![xy(&p.outline)];
                    rings.extend(p.holes.iter().filter(|h| h.len() >= 3).map(|h| xy(h)));
                    let poly = Polygon::new(rings);
                    if poly.area() > 0.0 {
                        out.push(CopperPoly {
                            net: p.net.clone(),
                            kind: p.kind.clone(),
                            poly,
                        });
                    }
                }
            }
            polys.push(out);
        }

        let pads = g
            .pads
            .into_iter()
            .map(|p| Pad {
                x: p.x - x0,
                y: -(p.y - y0),
                reference: p.reference,
                num: p.num,
                net: p.net,
                layers: p.layers,
            })
            .collect();
        let vias = g
            .vias
            .into_iter()
            .map(|v| Via {
                x: v.x - x0,
                y: -(v.y - y0),
                r: v.dia / 2.0,
                drill: v.drill.unwrap_or(0.0),
                top: v.top,
                bottom: v.bottom,
                net: v.net,
                kind: v.kind,
            })
            .collect();
        let footprints = g
            .footprints
            .into_iter()
            .map(|f| Footprint {
                reference: f.reference,
                value: f.value,
            })
            .collect();

        Ok(Board {
            board: g.board,
            kicad: g.kicad,
            layers: g.layers,
            x0,
            y0,
            w,
            h,
            cu_thickness: cu_t,
            z,
            diel,
            z_bot: zc,
            stackup_source,
            outline,
            polys,
            pads,
            vias,
            footprints,
        })
    }

    pub fn layer_index(&self, name: &str) -> Option<usize> {
        self.layers.iter().position(|l| l == name)
    }
}
