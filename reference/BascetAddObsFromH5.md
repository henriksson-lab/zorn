# Append Bascet HDF5/AnnData obs metadata to a Seurat object

Cells present in the Seurat object but missing from the HDF5 obs table
receive NA values. Cells present only in the HDF5 obs table are dropped.
Existing metadata columns with the same names are replaced.

## Usage

``` r
BascetAddObsFromH5(
  adata,
  fname,
  columns = NULL,
  prefix = NULL,
  subsetCommon = FALSE
)
```

## Arguments

- adata:

  A Seurat object

- fname:

  Path to a Bascet HDF5/AnnData file

- columns:

  Obs columns to add. Default: all columns except \_index

- prefix:

  Optional prefix added to column names before appending

- subsetCommon:

  If TRUE, subset the Seurat object to only cells present in both adata
  and obs

## Value

The Seurat object with added metadata columns
