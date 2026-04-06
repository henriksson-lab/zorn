# Aggregate data from AMRfinder This is a thin wrapper around BascetAggregateMap

Aggregate data from AMRfinder This is a thin wrapper around
BascetAggregateMap

## Usage

``` r
BascetAggregateAMRfinder(
  bascetRoot,
  inputName = "AMRfinder",
  includeCells = NULL,
  getColumn = "Element.symbol",
  ...
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- includeCells:

  Character vector of cell names to include, or NULL for all cells

- getColumn:

  Column name from AMRfinder output to use as feature (string)

- ...:

  Additional arguments passed to
  [`BascetAggregateMap`](https://henriksson-lab.github.io/zorn/reference/BascetAggregateMap.md)

## Value

Aggregated data

## See also

[`BascetAggregateMap`](https://henriksson-lab.github.io/zorn/reference/BascetAggregateMap.md)
