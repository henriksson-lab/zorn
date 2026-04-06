# Aggregate data from QUAST This is a thin wrapper around BascetAggregateMap

Aggregate data from QUAST This is a thin wrapper around
BascetAggregateMap

## Usage

``` r
BascetAggregateQUAST(bascetRoot, inputName = "quast", includeCells = NULL, ...)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- includeCells:

  Character vector of cell names to include, or NULL for all cells

- ...:

  Additional arguments passed to
  [`BascetAggregateMap`](https://henriksson-lab.github.io/zorn/reference/BascetAggregateMap.md)

## Value

Aggregated data

## See also

[`BascetAggregateMap`](https://henriksson-lab.github.io/zorn/reference/BascetAggregateMap.md)
