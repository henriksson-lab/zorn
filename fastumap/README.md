# fastumap

GPU-accelerated UMAP in Rust using CUDA. Computes 2D embeddings from high-dimensional data using cosine similarity.

Under testing

## Requirements

- CUDA toolkit (tested with 12.8)
- NVIDIA GPU
- Rust toolchain

## Build

```bash
cargo build --release
```

## Usage

```bash
./target/release/fastumap --input data.csv --output umap.csv
```

### Input format

CSV with no header:
- Column 1: cell/sample name (string)
- Column 2: depth (ignored)
- Columns 3+: numeric features

### Output format

CSV with header: `cell_name,umap_x,umap_y`

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--input` | required | Input CSV path |
| `--output` | required | Output CSV path |
| `--n-neighbors` | 15 | Number of nearest neighbors |
| `--min-dist` | 0.1 | Minimum distance in embedding |
| `--n-epochs` | 500 | Number of optimization epochs |
| `--learning-rate` | 1.0 | SGD learning rate |
| `--seed` | 42 | Random seed |
| `--gpu-device` | 0 | CUDA device index |

## How it works

1. **k-NN graph** — Data is uploaded to GPU, L2-normalized, then cosine similarity is computed via cuBLAS `sgemm` (matrix multiply). Top-k neighbors are extracted on CPU.
2. **Fuzzy simplicial set** — Membership strengths are computed per the UMAP algorithm (sigma binary search, symmetrization).
3. **Spectral initialization** — Initial 2D layout from the graph Laplacian via power iteration.
4. **SGD optimization** — Attractive/repulsive forces refine the layout over `n_epochs` iterations.

## Visualization

```r
source("cuda_umap.R")
```
