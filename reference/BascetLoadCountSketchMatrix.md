# Load count sketch matrix as Seurat object

Load count sketch matrix as Seurat object

## Usage

``` r
BascetLoadCountSketchMatrix(bascetRoot, inputName = "countsketch_mat.feather")
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of countsketch matrix file. Feather and legacy CSV files are
  supported.

## Value

A Seurat object holding the sketch as a reduction
