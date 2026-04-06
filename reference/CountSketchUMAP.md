# Run UMAP on a count sketch reduction

Run UMAP on a count sketch reduction

## Usage

``` r
CountSketchUMAP(adata, reduction = "kmersketch", metric = "cosine", ...)
```

## Arguments

- adata:

  A Seurat object with a count sketch reduction

- reduction:

  Name of the reduction to use

## Value

Seurat object with UMAP computed
