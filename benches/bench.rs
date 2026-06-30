use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_sc_pearson_residuals::{CountMatrix, Entry, PearsonParams, pearson_residuals};

fn make_matrix(n_cells: usize, n_genes: usize, density: f64, seed: u64) -> CountMatrix {
    // Simple LCG for reproducible sparse fill without external deps.
    let mut rng = seed;
    let lcg = |s: &mut u64| -> f64 {
        *s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (*s >> 33) as f64 / (u32::MAX as f64)
    };
    let mut entries = Vec::new();
    for ci in 0..n_cells {
        for gj in 0..n_genes {
            if lcg(&mut rng) < density {
                let count = (lcg(&mut rng) * 49.0) as i64 + 1;
                entries.push(Entry {
                    gene: gj as u32,
                    cell: ci as u32,
                    count,
                });
            }
        }
    }
    CountMatrix {
        n_genes,
        n_cells,
        entries,
    }
}

fn bench_compute(c: &mut Criterion) {
    // Representative single-cell dataset: 2000 cells × 3000 genes, ~10% density.
    let m = make_matrix(2000, 3000, 0.10, 42);
    let params = PearsonParams {
        theta: 100.0,
        clip: None,
    };

    c.bench_function("pearson_residuals_2000x3000", |b| {
        b.iter(|| {
            let res = pearson_residuals(black_box(&m), black_box(&params)).unwrap();
            black_box(res);
        });
    });
}

criterion_group!(benches, bench_compute);
criterion_main!(benches);
