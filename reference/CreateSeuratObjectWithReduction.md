# Create a seurat object from e.g. count sketch reduction

Create a seurat object from e.g. count sketch reduction

## Usage

``` r
CreateSeuratObjectWithReduction(Q, reductionName = "kmersketch", assay = "RNA")
```

## Arguments

- Q:

  Countsketch matrix

- reductionName:

  Name of reduction to store Q into

- assay:

  Name of assay to put into reduction

## Value

Seurat object holding Q as a reduction
