/// Compat tests against frozen scanpy 1.12.1 golden outputs.
///
/// Goldens were generated from real scanpy (not from this crate) and are stored
/// bit-exact in tests/golden/goldens.json. No Python is run at test time.
use std::fs;

use rsomics_sc_pearson_residuals::{CountMatrix, Entry, PearsonParams, pearson_residuals};
use serde::Deserialize;

#[derive(Deserialize)]
struct Goldens {
    scanpy_version: String,
    cases: Vec<Case>,
}

#[derive(Deserialize)]
struct Case {
    name: String,
    n_cells: usize,
    n_genes: usize,
    counts_flat_i64: Vec<i64>,
    theta: f64,
    clip: f64,
    residuals_flat_hex: Vec<String>,
}

fn hex_to_f64(h: &str) -> f64 {
    let bytes: [u8; 8] = hex::decode(h).unwrap().try_into().unwrap();
    f64::from_be_bytes(bytes)
}

fn matrix_from_flat(n_cells: usize, n_genes: usize, counts: &[i64]) -> CountMatrix {
    let mut entries = Vec::new();
    for ci in 0..n_cells {
        for gj in 0..n_genes {
            let v = counts[ci * n_genes + gj];
            if v != 0 {
                entries.push(Entry {
                    gene: gj as u32,
                    cell: ci as u32,
                    count: v,
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

#[test]
fn compat_all_cases() {
    let json_path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden/goldens.json");
    let raw = fs::read_to_string(json_path).expect("goldens.json missing");
    let goldens: Goldens = serde_json::from_str(&raw).expect("parse goldens.json");

    // Record provenance in output so CI logs capture it.
    eprintln!("goldens from scanpy {}", goldens.scanpy_version);

    for case in &goldens.cases {
        let m = matrix_from_flat(case.n_cells, case.n_genes, &case.counts_flat_i64);
        let params = PearsonParams {
            theta: case.theta,
            clip: Some(case.clip),
        };
        let got =
            pearson_residuals(&m, &params).unwrap_or_else(|e| panic!("case {}: {e}", case.name));

        assert_eq!(
            got.len(),
            case.residuals_flat_hex.len(),
            "case {}: output length mismatch",
            case.name
        );

        let mut max_err = 0.0f64;
        let mut nan_mismatches = 0usize;

        for (idx, (g, hex)) in got.iter().zip(case.residuals_flat_hex.iter()).enumerate() {
            let exp = hex_to_f64(hex);
            if exp.is_nan() {
                if !g.is_nan() {
                    nan_mismatches += 1;
                    eprintln!("case {}: idx={idx} expected NaN got {g}", case.name);
                }
            } else {
                let err = (g - exp).abs();
                if err > max_err {
                    max_err = err;
                }
            }
        }

        assert_eq!(
            nan_mismatches, 0,
            "case {}: {nan_mismatches} NaN position mismatches",
            case.name
        );
        assert!(
            max_err < 1e-13,
            "case {}: max |err| = {max_err:.2e} exceeds 1e-13",
            case.name
        );
        eprintln!("case {}: OK, max_err={max_err:.2e}", case.name);
    }
}
