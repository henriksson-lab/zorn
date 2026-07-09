# Bascet count matrix

Sparse count matrix container used by Zorn. Rows are observations/cells
and columns are features.

## Usage

``` r
# S4 method for class 'BascetCountMatrix'
show(object)

# S4 method for class 'BascetCountMatrix,ANY,ANY,ANY'
x[i, j, ..., drop = FALSE]
```

## Arguments

- object:

  A BascetCountMatrix object.

- x:

  A BascetCountMatrix object

- i:

  Row (cell) indices: numeric, logical, or character vector

- j:

  Column (feature) indices: numeric, logical, or character vector

- ...:

  Additional arguments ignored by this method

- drop:

  Ignored (kept for S4 signature compatibility)

## Value

A BascetCountMatrix

## Methods (by generic)

- `show(BascetCountMatrix)`: Print a compact summary of a Bascet count
  matrix.

- `x[i`: Subset rows and/or columns.

## Slots

- `X`:

  Sparse count matrix.

- `obs`:

  Data frame with observation metadata, one row per matrix row.
