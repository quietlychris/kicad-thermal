//! Results are bit-for-bit identical for any thread count (fixed chunking, ordered reductions).

mod common;
use common::*;
use kicad_thermal::{Grid, Precond};

#[test]
fn same_bits_on_1_and_4_threads() {
    let (bd, r) = load_case("b2_compact");
    let case = r.case().clone();
    for pc in [Precond::Jacobi, Precond::Amg] {
        let run_with = |threads: usize| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            pool.install(|| {
                let grid = Grid::new(&bd, case["cell"].as_f64().unwrap()).unwrap();
                let run = run_case(&grid, &case, pc);
                (run.th.t, run.dc.unwrap().v)
            })
        };
        let (t1, v1) = run_with(1);
        let (t4, v4) = run_with(4);
        assert!(
            t1.iter().zip(&t4).all(|(a, b)| a.to_bits() == b.to_bits()),
            "{pc:?}: T differs between 1 and 4 threads"
        );
        assert!(
            v1.iter().zip(&v4).all(|(a, b)| a.to_bits() == b.to_bits()),
            "{pc:?}: V differs between 1 and 4 threads"
        );
    }
}
