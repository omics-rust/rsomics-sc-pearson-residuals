use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};

pub struct CountMatrix {
    pub n_genes: usize,
    pub n_cells: usize,
    /// Nonzero entries from the MTX file. 10x stores genes on rows, cells on cols.
    pub entries: Vec<Entry>,
}

#[derive(Clone, Copy)]
pub struct Entry {
    pub gene: u32,
    pub cell: u32,
    pub count: i64,
}

pub struct PearsonParams {
    /// Negative-binomial overdispersion shared across genes. Must be > 0.
    pub theta: f64,
    /// Symmetric clip bound. `None` = sqrt(n_cells) per scanpy default.
    pub clip: Option<f64>,
}

pub fn open_mtx(dir: &Path) -> Result<Box<dyn Read>> {
    for name in ["matrix.mtx.gz", "matrix.mtx"] {
        let path = dir.join(name);
        if path.exists() {
            return open_maybe_gz(&path);
        }
    }
    Err(RsomicsError::InvalidInput(format!(
        "no matrix.mtx or matrix.mtx.gz in {}",
        dir.display()
    )))
}

fn open_maybe_gz(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

/// Parse a MatrixMarket coordinate file. 10x stores genes on rows, cells on cols.
pub fn parse_mtx(reader: impl Read) -> Result<CountMatrix> {
    let mut reader = BufReader::new(reader);
    let mut line = String::new();

    reader.read_line(&mut line).map_err(RsomicsError::Io)?;
    let banner = line.trim();
    if !banner.starts_with("%%MatrixMarket") {
        return Err(RsomicsError::InvalidInput(
            "missing %%MatrixMarket banner".into(),
        ));
    }
    let pattern = banner.contains("pattern");

    let (n_genes, n_cells, nnz) = loop {
        line.clear();
        let n = reader.read_line(&mut line).map_err(RsomicsError::Io)?;
        if n == 0 {
            return Err(RsomicsError::InvalidInput("truncated MTX header".into()));
        }
        let t = line.trim();
        if t.is_empty() || t.starts_with('%') {
            continue;
        }
        let mut it = t.split_whitespace();
        let rows = parse_usize(it.next())?;
        let cols = parse_usize(it.next())?;
        let nnz = parse_usize(it.next())?;
        break (rows, cols, nnz);
    };

    let mut entries = Vec::with_capacity(nnz);
    for raw in reader.lines() {
        let raw = raw.map_err(RsomicsError::Io)?;
        let t = raw.trim();
        if t.is_empty() {
            continue;
        }
        let mut it = t.split_whitespace();
        let gene = parse_usize(it.next())?;
        let cell = parse_usize(it.next())?;
        let count: i64 = if pattern {
            1
        } else {
            it.next()
                .ok_or_else(|| RsomicsError::InvalidInput("MTX entry missing value".into()))?
                .parse::<i64>()?
        };
        if gene == 0 || gene > n_genes || cell == 0 || cell > n_cells {
            return Err(RsomicsError::InvalidInput(format!(
                "MTX index out of bounds: ({gene}, {cell})"
            )));
        }
        entries.push(Entry {
            gene: (gene - 1) as u32,
            cell: (cell - 1) as u32,
            count,
        });
    }
    if entries.len() != nnz {
        return Err(RsomicsError::InvalidInput(format!(
            "MTX declared {nnz} entries, found {}",
            entries.len()
        )));
    }

    Ok(CountMatrix {
        n_genes,
        n_cells,
        entries,
    })
}

/// Compute analytic Pearson residuals, returning a dense `n_cells × n_genes`
/// matrix (row-major, cell-major outer index) matching scanpy's formula:
///
///   mu_ij = row_sum_i * col_sum_j / total
///   r_ij  = (x_ij - mu_ij) / sqrt(mu_ij + mu_ij^2 / theta)
///   clipped to [-clip, +clip]  where clip defaults to sqrt(n_cells)
///
/// Integer row/col/total sums are exact; residuals then use f64 arithmetic.
/// When mu_ij = 0 (zero-count gene or zero-count cell), the residual is NaN,
/// matching scanpy's 0/0 behaviour. Callers that need to exclude such genes
/// should filter them before calling.
pub fn pearson_residuals(m: &CountMatrix, params: &PearsonParams) -> Result<Vec<f64>> {
    if params.theta <= 0.0 {
        return Err(RsomicsError::InvalidInput(
            "Pearson residuals require theta > 0".into(),
        ));
    }

    let n_cells = m.n_cells;
    let n_genes = m.n_genes;
    let clip = params.clip.unwrap_or_else(|| (n_cells as f64).sqrt());

    // Integer row (cell) and column (gene) sums for exact accumulation.
    let mut row_sum = vec![0i64; n_cells];
    let mut col_sum = vec![0i64; n_genes];
    for e in &m.entries {
        row_sum[e.cell as usize] += e.count;
        col_sum[e.gene as usize] += e.count;
    }
    let total: i64 = col_sum.iter().sum();

    // Precompute float versions used in the inner loops.
    let total_f = total as f64;
    let row_f: Vec<f64> = row_sum.iter().map(|&s| s as f64).collect();
    let col_f: Vec<f64> = col_sum.iter().map(|&s| s as f64).collect();

    // Dense residual buffer: row-major, cell-major (n_cells × n_genes).
    // Initialise from the implicit-zero contribution: x_ij=0,
    //   r_ij = (0 - mu_ij) / denom_ij = -mu_ij / denom_ij
    let mut dense = vec![0.0f64; n_cells * n_genes];
    // A gene-less matrix has an empty result; par_chunks_mut(0) would panic,
    // and scanpy returns the defined (n_cells, 0) array here.
    if n_genes == 0 {
        return Ok(dense);
    }
    dense
        .par_chunks_mut(n_genes)
        .enumerate()
        .for_each(|(ci, row)| {
            let rs = row_f[ci];
            for gj in 0..n_genes {
                let mu = rs * col_f[gj] / total_f;
                let denom = (mu + mu * mu / params.theta).sqrt();
                // 0/0 when mu=0: matches scanpy's NaN propagation.
                row[gj] = (-mu / denom).clamp(-clip, clip);
            }
        });

    // Overwrite the nonzero positions.
    for e in &m.entries {
        let ci = e.cell as usize;
        let gj = e.gene as usize;
        let mu = row_f[ci] * col_f[gj] / total_f;
        let denom = (mu + mu * mu / params.theta).sqrt();
        let r = (e.count as f64 - mu) / denom;
        dense[ci * n_genes + gj] = r.clamp(-clip, clip);
    }

    Ok(dense)
}

/// Write the dense residual matrix as MatrixMarket array real general layout
/// (column-major: gene 0 for all cells, gene 1 for all cells, …), matching
/// scipy's dense MM writer that scanpy users pipe into.
pub fn write_dense(n_cells: usize, n_genes: usize, dense: &[f64], out: impl Write) -> Result<()> {
    // MM array is column-major. Our dense is row-major (cell × gene).
    // Column j of MM = gene j across all cells = column j of dense.
    let mut w = BufWriter::with_capacity(1 << 20, out);
    w.write_all(b"%%MatrixMarket matrix array real general\n")
        .map_err(RsomicsError::Io)?;
    writeln!(w, "{n_cells} {n_genes}").map_err(RsomicsError::Io)?;

    let mut fmt = ryu::Buffer::new();
    let mut buf: Vec<u8> = Vec::with_capacity(1 << 16);
    // Column-major order: iterate genes first.
    for gj in 0..n_genes {
        for ci in 0..n_cells {
            buf.extend_from_slice(fmt.format(dense[ci * n_genes + gj]).as_bytes());
            buf.push(b'\n');
            if buf.len() >= 1 << 15 {
                w.write_all(&buf).map_err(RsomicsError::Io)?;
                buf.clear();
            }
        }
    }
    w.write_all(&buf).map_err(RsomicsError::Io)?;
    w.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

pub fn open_output(path: &str) -> Result<Box<dyn Write>> {
    if path == "-" {
        Ok(Box::new(std::io::stdout().lock()))
    } else {
        Ok(Box::new(
            File::create(PathBuf::from(path)).map_err(RsomicsError::Io)?,
        ))
    }
}

fn parse_usize(tok: Option<&str>) -> Result<usize> {
    tok.ok_or_else(|| RsomicsError::InvalidInput("MTX header missing a dimension".into()))?
        .parse::<usize>()
        .map_err(Into::into)
}

/// End-to-end: read 10x MTX, compute Pearson residuals, write dense MM.
pub fn run(dir: &Path, params: &PearsonParams, out: impl Write) -> Result<(usize, usize)> {
    let m = parse_mtx(open_mtx(dir)?)?;
    let shape = (m.n_genes, m.n_cells);
    let dense = pearson_residuals(&m, params)?;
    write_dense(m.n_cells, m.n_genes, &dense, out)?;
    Ok(shape)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix_from_dense(counts: &[&[i64]]) -> CountMatrix {
        // counts[cell][gene]
        let n_cells = counts.len();
        let n_genes = counts[0].len();
        let mut entries = Vec::new();
        for (ci, row) in counts.iter().enumerate() {
            for (gj, &v) in row.iter().enumerate() {
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
    fn matches_scanpy_small_dense() {
        // 4x5 from goldens (rng seed 42):
        // [[1,15,13,8,8],[17,1,13,4,1],[10,19,14,15,14],[15,10,2,16,9]]
        // scanpy residuals (row-major):
        // [[-2.0, 1.5546921, 1.1913620, -0.4477324, 0.3558201],
        //  [ 2.0,-2.0,       1.9985898, -1.2461241,-1.8961589],
        //  [-1.2238044, 0.7468415,-0.1825889,-0.0245697, 0.7808537],
        //  [ 1.1767089,-0.3966816,-2.0,       1.4642242, 0.2980423]]
        let counts: &[&[i64]] = &[
            &[1, 15, 13, 8, 8],
            &[17, 1, 13, 4, 1],
            &[10, 19, 14, 15, 14],
            &[15, 10, 2, 16, 9],
        ];
        let m = matrix_from_dense(counts);
        let params = PearsonParams {
            theta: 100.0,
            clip: Some(2.0),
        };
        let res = pearson_residuals(&m, &params).unwrap();

        let expected: &[f64] = &[
            -2.0,
            1.55469208,
            1.19136203,
            -0.44773237,
            0.35582007,
            2.0,
            -2.0,
            1.99858977,
            -1.24612412,
            -1.89615891,
            -1.22380441,
            0.74684149,
            -0.18258885,
            -0.02456969,
            0.78085371,
            1.17670888,
            -0.39668157,
            -2.0,
            1.46422416,
            0.29804226,
        ];
        for (got, &exp) in res.iter().zip(expected.iter()) {
            assert!((got - exp).abs() < 1e-6, "got {got} expected {exp}");
        }
    }

    #[test]
    fn clip_default_is_sqrt_ncells() {
        let counts: &[&[i64]] = &[
            &[1, 15, 13, 8, 8],
            &[17, 1, 13, 4, 1],
            &[10, 19, 14, 15, 14],
            &[15, 10, 2, 16, 9],
        ];
        let m = matrix_from_dense(counts);
        let params = PearsonParams {
            theta: 100.0,
            clip: None,
        };
        let res = pearson_residuals(&m, &params).unwrap();
        let clip = (4.0f64).sqrt();
        for &v in &res {
            assert!(v.abs() <= clip + 1e-12, "v={v} exceeds clip={clip}");
        }
    }

    #[test]
    fn zero_gene_column_produces_nan() {
        // Gene 1 has zero counts in all cells -> mu=0 -> NaN residual.
        let counts: &[&[i64]] = &[&[3, 0, 1], &[5, 0, 2], &[1, 0, 3]];
        let m = matrix_from_dense(counts);
        let params = PearsonParams {
            theta: 100.0,
            clip: None,
        };
        let res = pearson_residuals(&m, &params).unwrap();
        // cells × genes = 3×3, gene index 1 positions: res[0*3+1], res[1*3+1], res[2*3+1]
        for ci in 0..3 {
            assert!(res[ci * 3 + 1].is_nan(), "expected NaN for zero-gene col");
        }
        // Non-zero genes should be finite.
        for ci in 0..3 {
            for gj in [0usize, 2] {
                assert!(
                    res[ci * 3 + gj].is_finite(),
                    "expected finite for nonzero gene"
                );
            }
        }
    }

    #[test]
    fn theta_gt_zero_required() {
        let m = CountMatrix {
            n_genes: 1,
            n_cells: 1,
            entries: vec![],
        };
        assert!(
            pearson_residuals(
                &m,
                &PearsonParams {
                    theta: 0.0,
                    clip: None
                }
            )
            .is_err()
        );
        assert!(
            pearson_residuals(
                &m,
                &PearsonParams {
                    theta: -1.0,
                    clip: None
                }
            )
            .is_err()
        );
    }

    #[test]
    fn zero_gene_matrix_returns_empty() {
        // scanpy returns a defined (n_cells, 0) array; we must not panic in
        // par_chunks_mut(0).
        let m = CountMatrix {
            n_genes: 0,
            n_cells: 5,
            entries: vec![],
        };
        let res = pearson_residuals(
            &m,
            &PearsonParams {
                theta: 100.0,
                clip: None,
            },
        )
        .unwrap();
        assert!(res.is_empty());
    }

    #[test]
    fn zero_by_zero_matrix_returns_empty() {
        let m = CountMatrix {
            n_genes: 0,
            n_cells: 0,
            entries: vec![],
        };
        let res = pearson_residuals(
            &m,
            &PearsonParams {
                theta: 100.0,
                clip: None,
            },
        )
        .unwrap();
        assert!(res.is_empty());
    }

    #[test]
    fn roundtrip_mtx_parse() {
        let counts: &[&[i64]] = &[&[3, 2, 0], &[0, 5, 1]];
        let m = matrix_from_dense(counts);
        let mut buf = Vec::new();
        // Write as coordinate MTX manually.
        buf.extend_from_slice(b"%%MatrixMarket matrix coordinate integer general\n");
        buf.extend_from_slice(b"2 3 4\n"); // n_genes=2, n_cells=3 in MTX (genes on rows)
        // entry: gene(row) cell(col) count — 1-based
        buf.extend_from_slice(b"1 1 3\n1 2 2\n2 2 5\n2 3 1\n");
        let parsed = parse_mtx(&buf[..]).unwrap();
        assert_eq!(parsed.n_genes, 2);
        assert_eq!(parsed.n_cells, 3);
        assert_eq!(parsed.entries.len(), 4);
        let _ = m; // suppressed
    }
}
