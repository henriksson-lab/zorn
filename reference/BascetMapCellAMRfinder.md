# Run AMRfinder on contigs of all cells. This is a thin wrapper around BascetMapCell

Run AMRfinder on contigs of all cells. This is a thin wrapper around
BascetMapCell

## Usage

``` r
BascetMapCellAMRfinder(
  bascetRoot,
  inputName = "contigs",
  outputName = "AMRfinder",
  db,
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

- db:

  Path to database

- ...:

  Additional arguments passed to
  [`BascetMapCell`](https://henriksson-lab.github.io/zorn/reference/BascetMapCell.md)

## Value

A job to be executed, or being executed, depending on runner settings

## See also

[`BascetMapCell`](https://henriksson-lab.github.io/zorn/reference/BascetMapCell.md)
