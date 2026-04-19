# Run Bakta on contigs of all cells. This is a thin wrapper around BascetMapCell

Run Bakta on contigs of all cells. This is a thin wrapper around
BascetMapCell

## Usage

``` r
BascetMapCellBakta(
  bascetRoot,
  inputName = "contigs",
  outputName = "bakta",
  db,
  args = list(),
  overwrite = FALSE,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
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

- args:

  Additional arguments passed to Bascet (list)

- overwrite:

  Whether to overwrite existing output (logical)

- runner:

  A runner object (e.g. LocalRunner, SlurmRunner)

- bascetInstance:

  A Bascet instance

## Value

A job to be executed, or being executed, depending on runner settings

## See also

[`BascetMapCell`](https://henriksson-lab.github.io/zorn/reference/BascetMapCell.md)
