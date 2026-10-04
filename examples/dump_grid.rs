//! Debug helper: write a grid's materials as raw little-endian f64 for comparison with Python.
//!     cargo run --release --example dump_grid -- geometry.json 0.25 outdir
use kicad_thermal::{Board, Grid};
use std::io::Write;

fn write_f64(path: &std::path::Path, v: &[f64]) {
    let mut f = std::fs::File::create(path).unwrap();
    for x in v {
        f.write_all(&x.to_le_bytes()).unwrap();
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let bd = Board::load(&a[1]).unwrap();
    let g = Grid::new(&bd, a[2].parse().unwrap()).unwrap();
    let out = std::path::Path::new(&a[3]);
    std::fs::create_dir_all(out).unwrap();
    write_f64(&out.join("kx.f64"), &g.kx);
    write_f64(&out.join("kz.f64"), &g.kz);
    write_f64(&out.join("via.f64"), &g.via_copper_area(None));
    for (l, f) in g.cu_frac.iter().enumerate() {
        write_f64(&out.join(format!("cufrac_{l}.f64")), f);
    }
    let inside: Vec<f64> = g.inside.iter().map(|&b| b as u8 as f64).collect();
    write_f64(&out.join("inside.f64"), &inside);
    println!("{} {} {}", g.nx, g.ny, g.nz);
}
