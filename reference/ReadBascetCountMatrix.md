# Read a count matrix as produced by Bascet (hdf5 format). This can be output from both BascetQueryFq and BascetCountChrom

Read a count matrix as produced by Bascet (hdf5 format). This can be
output from both BascetQueryFq and BascetCountChrom

## Usage

``` r
ReadBascetCountMatrix(bascetRoot, inputName, verbose = FALSE)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- verbose:

  Print additional information, primarily to help troubleshooting

## Value

Count matrix as sparseMatrix
