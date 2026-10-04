//! Shared helpers: fixture paths, reference data loading, case runner.
#![allow(dead_code)]

use flate2::read::GzDecoder;
use kicad_thermal::{solve_dc, solve_thermal, Board, DcResult, Grid, Precond, ThermalResult};
use serde_json::Value;
use std::io::Read;
use std::path::PathBuf;

pub fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests")
}

pub fn fixture(name: &str) -> PathBuf {
    root().join("fixtures").join(format!("{name}.json"))
}

pub struct Reference {
    pub dir: PathBuf,
    pub meta: Value,
}

impl Reference {
    pub fn load(case: &str) -> Reference {
        let dir = root().join("reference").join(case);
        let text = std::fs::read_to_string(dir.join("meta.json")).unwrap_or_else(|e| {
            panic!("missing reference for {case} ({e}); run scripts/gen_reference.py")
        });
        Reference {
            dir,
            meta: serde_json::from_str(&text).unwrap(),
        }
    }

    pub fn array_file(&self, file: &str) -> Vec<f64> {
        let mut bytes = Vec::new();
        GzDecoder::new(std::fs::File::open(self.dir.join(file)).unwrap())
            .read_to_end(&mut bytes)
            .unwrap();
        bytes
            .chunks_exact(8)
            .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
            .collect()
    }

    /// A grid-level array ("kx", "kz", "via_area", "inside", "cu_frac_<l>").
    pub fn grid_array(&self, key: &str) -> Vec<f64> {
        self.array_file(self.meta["arrays"][key].as_str().unwrap())
    }

    /// A result array ("T", "V", "J") for a Python preconditioner ("jacobi" or "amg").
    pub fn result_array(&self, precond: &str, key: &str) -> Option<Vec<f64>> {
        self.meta["results"][precond]["arrays"][key]
            .as_str()
            .map(|f| self.array_file(f))
    }

    pub fn result(&self, precond: &str) -> &Value {
        &self.meta["results"][precond]
    }

    pub fn case(&self) -> &Value {
        &self.meta["case"]
    }
}

pub struct Run {
    pub dc: Option<DcResult>,
    pub th: ThermalResult,
    pub report: String,
}

/// Run one case from cases.json with a given preconditioner.
pub fn run_case(grid: &Grid, case: &Value, precond: Precond) -> Run {
    let dc = case["dc"].as_object().map(|d| {
        let s = |v: &Value, i: usize| v[i].as_str().unwrap().to_string();
        let (f, t) = (&d["from"], &d["to"]);
        let (f0, f1, t0, t1) = (s(f, 0), s(f, 1), s(t, 0), s(t, 1));
        solve_dc(
            grid,
            d["net"].as_str().unwrap(),
            (&f0, &f1),
            (&t0, &t1),
            d["current"].as_f64().unwrap(),
            precond,
        )
        .unwrap()
    });
    let powers: Vec<(String, f64)> = case["powers"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| (p[0].as_str().unwrap().to_string(), p[1].as_f64().unwrap()))
        .collect();
    let th = solve_thermal(
        grid,
        &powers,
        case["h"].as_f64().unwrap(),
        case["tamb"].as_f64().unwrap(),
        dc.as_ref().map(|d| d.joule.as_slice()),
        precond,
    )
    .unwrap();
    let report = kicad_thermal::report::report(grid, &th, dc.as_ref()).unwrap();
    Run { dc, th, report }
}

pub fn load_case(name: &str) -> (Board, Reference) {
    (Board::load(fixture(name)).unwrap(), Reference::load(name))
}

pub fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "field sizes differ");
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

pub fn rel(a: f64, b: f64) -> f64 {
    if a == b {
        0.0
    } else {
        (a - b).abs() / b.abs().max(f64::MIN_POSITIVE)
    }
}
