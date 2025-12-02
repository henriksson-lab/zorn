# Read a count matrix as produced by CellSNP, but as shards

Read a count matrix as produced by CellSNP, but as shards

## Usage

``` r
ReadCellSNPmatrix(bascetRoot, inputName, listCells = NULL, verbose = FALSE)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- listCells:

  Optional list of cells to extract

- verbose:

  Show process status for debugging purposes

## Value

Count matrix as sparseMatrix
