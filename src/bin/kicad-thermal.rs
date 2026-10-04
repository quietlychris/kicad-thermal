//! kicad-thermal: command-line equivalent of `kicad-fdtd thermal`.

use anyhow::{bail, Context, Result};
use clap::Parser;
use kicad_thermal::{cg::SolveStats, report, solve_dc, solve_thermal, Board, Grid, Precond};
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Steady-state thermal and DC IR drop on a board's exact copper.
///
/// Input is the geometry JSON written by `kicad-fdtd export board.kicad_pcb` (that step
/// needs KiCad's Python). Thermal: finite-volume heat conduction on a whole-board grid,
/// copper sampled 3 x 3 per cell, FR4 0.8 W/mK in plane and 0.3 through, via barrels,
/// component power under the pads, convection h to T_amb on both faces, edges adiabatic,
/// no radiation (a conservative estimate). DC: sheet and via conduction on one net,
/// current in at one pad and out at another; the Joule heat is added to the thermal sources.
#[derive(Parser, Debug)]
#[command(name = "kicad-thermal", version)]
struct Args {
    /// geometry JSON from `kicad-fdtd export`
    geometry: PathBuf,
    /// component power, REF=W ...
    #[arg(long, num_args = 0.., value_name = "REF=W")]
    power: Vec<String>,
    /// convection W/m2K on both faces
    #[arg(long, default_value_t = 10.0)]
    h: f64,
    /// ambient temperature, C
    #[arg(long, default_value_t = 25.0)]
    tamb: f64,
    /// in-plane cell size, mm
    #[arg(long, default_value_t = 0.25)]
    cell: f64,
    /// net for the DC solve
    #[arg(long)]
    net: Option<String>,
    /// DC source pad, REF.PAD
    #[arg(long = "from")]
    src: Option<String>,
    /// DC sink pad, REF.PAD
    #[arg(long = "to")]
    sink: Option<String>,
    /// DC current, A
    #[arg(long, default_value_t = 0.0)]
    current: f64,
    /// CG preconditioner
    #[arg(long, value_enum, default_value_t = Precond::Auto)]
    precond: Precond,
    /// accepted for compatibility with kicad-fdtd (this port always runs on the CPU)
    #[arg(long)]
    cpu: bool,
    /// output directory for the result JSON
    #[arg(long, default_value = "out")]
    out: PathBuf,
    /// print solver statistics to stderr
    #[arg(short, long)]
    verbose: bool,
}

fn split_pad(s: &str) -> Result<(&str, &str)> {
    s.split_once('.')
        .with_context(|| format!("pad '{s}' must be REF.PAD"))
}

fn log_stats(what: &str, s: &SolveStats) {
    let mut line = format!(
        "{what}: {} unknowns, {} preconditioner, {} iterations, setup {:.2} s, solve {:.2} s",
        s.n, s.precond, s.iterations, s.setup_s, s.solve_s
    );
    if let Some(a) = &s.amg {
        line += &format!(
            ", levels {:?}, operator complexity {:.2}",
            a.sizes,
            a.operator_complexity()
        );
    }
    if let Some(f) = &s.fallback {
        line += &format!(" ({f})");
    }
    eprintln!("{line}");
}

fn main() -> Result<()> {
    let a = Args::parse();
    if a.geometry.extension().and_then(|e| e.to_str()) == Some("kicad_pcb") {
        bail!("{} is a board file; export it first with `kicad-fdtd export {}` and pass the geometry JSON", a.geometry.display(), a.geometry.display());
    }
    let bd = Board::load(&a.geometry)?;
    let t = Instant::now();
    let g = Grid::new(&bd, a.cell)?;
    println!(
        "grid {} x {} x {} cells ({:.1} s)",
        g.nx,
        g.ny,
        g.nz,
        t.elapsed().as_secs_f64()
    );
    let dc = match &a.net {
        Some(net) => {
            let (Some(src), Some(sink)) = (&a.src, &a.sink) else {
                bail!("--net needs --from REF.PAD --to REF.PAD --current A");
            };
            if a.current == 0.0 {
                bail!("--net needs --from REF.PAD --to REF.PAD --current A");
            }
            let dc = solve_dc(
                &g,
                net,
                split_pad(src)?,
                split_pad(sink)?,
                a.current,
                a.precond,
            )?;
            if a.verbose {
                log_stats("DC", &dc.stats);
            }
            Some(dc)
        }
        None => None,
    };
    let mut powers = Vec::new();
    for s in &a.power {
        let (r, w) = s
            .split_once('=')
            .with_context(|| format!("--power entry '{s}' must be REF=W"))?;
        powers.push((
            r.to_string(),
            w.parse::<f64>()
                .with_context(|| format!("bad watts in '{s}'"))?,
        ));
    }
    let th = solve_thermal(
        &g,
        &powers,
        a.h,
        a.tamb,
        dc.as_ref().map(|d| d.joule.as_slice()),
        a.precond,
    )?;
    if a.verbose {
        log_stats("thermal", &th.stats);
    }
    println!("{}", report::report(&g, &th, dc.as_ref())?);
    std::fs::create_dir_all(&a.out)?;
    let stem = Path::new(&a.geometry)
        .file_name()
        .and_then(|f| f.to_str())
        .unwrap_or("board");
    let stem = stem.split('_').next().unwrap_or(stem);
    let path = a.out.join(format!("thermal_{stem}.json"));
    std::fs::write(
        &path,
        serde_json::to_string(&report::result_json(&g, &th, dc.as_ref())?)?,
    )?;
    println!("saved {}", path.display());
    Ok(())
}
