# rsomics-sc-pearson-residuals

Analytic Pearson residuals normalization of a single-cell count matrix.

Ports `scanpy.experimental.pp.normalize_pearson_residuals` (Lause et al., 2021). Reads a 10x
MatrixMarket count matrix (genes × cells), applies the negative-binomial-based residuals
formula, and writes a dense residual matrix. Distinct from `rsomics-sc-normalize`, which does
library-size + log1p normalization.

## Install

```
cargo install rsomics-sc-pearson-residuals
```

## Usage

```
rsomics-sc-pearson-residuals <10x-mtx-dir> [OPTIONS]

OPTIONS:
  -o, --output <path>      Output dense MM path; '-' for stdout (default: -)
      --theta <float>      NB overdispersion, shared across genes (default: 100.0)
      --clip <float>       Symmetric clip bound; default = sqrt(n_cells)
  -t, --threads <N>        Rayon thread count (default: all cores)
  -q, --quiet              Suppress progress message
```

The output is the residual matrix in MatrixMarket array format (the round-trippable
single-cell convention); `--json` is not supported and exits with an error rather than
silently producing the matrix.

The input directory must contain `matrix.mtx` or `matrix.mtx.gz` (genes on rows, cells on
columns — standard 10x layout). Output is MatrixMarket `array real general` (column-major),
cells × genes.

## Formula

For cell *i*, gene *j*:

```
mu_ij    = row_sum_i × col_sum_j / total
r_ij     = (x_ij − mu_ij) / sqrt(mu_ij + mu_ij² / theta)
clipped to [−clip, +clip]
```

Row/column/total sums accumulate in `i64` (exact). Division and sqrt use `f64`. When
`mu_ij = 0` (zero-count gene or zero-count cell), the residual is `NaN`, matching scanpy's
behaviour. Filter zero-count genes before calling if this is undesirable.

## Accuracy

Value-exact vs scanpy 1.12.1 across all test cases: **bit-exact (max |err| = 0)**. The integer
accumulation path eliminates floating-point rounding in the sums, matching scanpy's
`np.sum(x, axis=0)` integer path exactly when counts are integers.

## Origin

This crate is an independent Rust reimplementation based on:

- Lause, J., Berens, P. & Kobak, D. (2021). Analytic Pearson residuals for normalization of
  single-cell RNA-seq UMI data. *Genome Biology* **22**, 258.
  DOI: [10.1186/s13059-021-02451-7](https://doi.org/10.1186/s13059-021-02451-7)
- The scanpy BSD-3-Clause source (`scanpy/experimental/pp/_normalization.py`, version 1.12.1)
  was read for formula verification (BSD allows this).
- Golden test fixtures were generated from real scanpy 1.12.1 and are frozen in
  `tests/golden/goldens.json`. No Python runs at test time.

License: MIT OR Apache-2.0. Upstream: [scverse/scanpy](https://github.com/scverse/scanpy)
(BSD-3-Clause).
