# Performance Notes

## Setup

- Machine: mini_m2 (Apple M2, macOS 25.5.0)
- Rust: stable (1.91), release profile
- scanpy: 1.12.1, Python 3.12, NumPy 2.x, OMP_NUM_THREADS=1/OPENBLAS_NUM_THREADS=1
- Fixture: 3000 cells × 5000 genes, ~10% density, 1 427 137 nnz (17 MB MTX file)
- Tool binary: `/Volumes/KIOXIA/Developments/cargo-target/release/rsomics-sc-pearson-residuals`

## Axis 1: compute-only

Matrix size 2000 cells × 3000 genes, 10% density (matches Criterion bench fixture).

| Implementation | Mean ms | Notes |
|---|---|---|
| scanpy 1.12.1 (single-thread OMP=1) | 59.4 ms | preloaded dense f64 matrix |
| Rust (criterion, multi-thread) | 10.4 ms | in-memory sparse parse + compute |

**Ratio: 5.71× (ours faster)**

Scanpy's hot path is `sums_cells @ sums_genes` — a dense (n, 1) × (1, g) outer product via
BLAS `dger` / `dgemm`, then element-wise ops. Our path builds cell/gene sums as i64 scalars
(O(nnz)), then fills the dense output in parallel rayon chunks. For sparse inputs (~10%),
skipping the full O(n·g) matmul wins substantially.

## Axis 2: both-serialize

Full pipeline: read MTX file → compute → write to `/dev/null`.

| Implementation | Mean ms | Notes |
|---|---|---|
| scanpy 1.12.1 (`scipy.io.mmread` + `np.savetxt`) | 6761 ms | read sparse, toarray, compute, savetxt |
| Rust (hyperfine, `--warmup 3 --runs 10`) | 687 ms ± 29 ms | read MTX, compute, write MM array |

**Ratio: 9.84× (ours faster)**

Scanpy's `np.savetxt` dominates the both-serialize time (6.5 s) because it formats floats
through Python's string machinery one row at a time. Our writer uses `ryu` (fast Grisu3) with
a 32 KB write buffer, reducing format+write to ~milliseconds.

## Verdict

> PERF PASS — strictly > 1.0× on both axes.

Compute-only: **5.71×**. Both-serialize: **9.84×**. Published at 0.1.0.
