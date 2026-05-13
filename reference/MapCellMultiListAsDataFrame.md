# Convenience function; alternative is to somehow implement as.data.frame

This one puts cellID as a new column

## Usage

``` r
MapCellMultiListAsDataFrame(mylist)
```

## Arguments

- mylist:

  Named list of data.frames, keyed by cell ID

## Value

A data.frame containing the row-bound non-empty entries of `mylist`,
with a `cellID` column identifying the source list element.

## Details

FUTURE possible make this the new MapListAsDataFrame
