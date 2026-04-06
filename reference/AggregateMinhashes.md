# Aggregate frequency of minhashes across cells

Aggregate frequency of minhashes across cells

## Usage

``` r
AggregateMinhashes(bascetRoot, inputName = "minhash", includeCells = NULL, ...)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard (Container of minhashes)

- includeCells:

  Character vector of cell names to include, or NULL for all cells

- ...:

  Additional arguments passed to
  [`BascetAggregateMap`](https://henriksson-lab.github.io/zorn/reference/BascetAggregateMap.md)

## Value

Data.frame of KMERs and frequencies

## See also

[`BascetAggregateMap`](https://henriksson-lab.github.io/zorn/reference/BascetAggregateMap.md)
