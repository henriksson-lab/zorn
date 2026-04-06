# Run SKESA on reads of all cells. This is a thin wrapper around BascetMapCell

Run SKESA on reads of all cells. This is a thin wrapper around
BascetMapCell

## Usage

``` r
BascetMapCellSKESA(
  bascetRoot,
  inputName = "filtered",
  outputName = "contigs",
  ...
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- outputName:

  Name of output shard

- ...:

  Additional arguments passed to
  [`BascetMapCell`](https://henriksson-lab.github.io/zorn/reference/BascetMapCell.md)

## Value

A job to be executed, or being executed, depending on runner settings

## See also

[`BascetMapCell`](https://henriksson-lab.github.io/zorn/reference/BascetMapCell.md)
