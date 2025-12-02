# List files for a cell in a Bascet

This can be made faster by, e.g., once and for all reading the location
of all objects in the file

## Usage

``` r
BascetListFilesForCell(
  bascetFile,
  cellID,
  bascetInstance = GetDefaultBascetInstance(),
  superVerbose = FALSE
)
```

## Arguments

- bascetFile:

  Bascet file object

- cellID:

  Name of the cell

- bascetInstance:

  A Bascet instance

## Value

A data.frame with list of all the files
