# Run UMAP on a count sketch reduction

Run UMAP on a count sketch reduction

## Usage

``` r
CountSketchUMAP(
  adata,
  reduction = "kmersketch",
  outReduction = "kmersketch_umap",
  metric = "cosine",
  n_neighbors = 30,
  doFast = FALSE,
  seed = 42
)
```

## Arguments

- adata:

  A Seurat object with a count sketch reduction

- reduction:

  Name of the reduction to use

- outReduction:

  Name of the output UMAP reduction

- metric:

  Distance metric passed to
  [`uwot::umap`](https://jlmelville.github.io/uwot/reference/umap.html)

- n_neighbors:

  Number of nearest neighbors passed to
  [`uwot::umap`](https://jlmelville.github.io/uwot/reference/umap.html)

- doFast:

  Use multi-threaded stochastic gradient updates

- seed:

  Random seed passed to
  [`uwot::umap`](https://jlmelville.github.io/uwot/reference/umap.html)

## Value

Seurat object with UMAP computed
