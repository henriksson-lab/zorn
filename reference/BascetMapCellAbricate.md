# Run Abricate on contigs of all cells. This is a thin wrapper around BascetMapCell

Run Abricate on contigs of all cells. This is a thin wrapper around
BascetMapCell

## Usage

``` r
BascetMapCellAbricate(
  bascetRoot,
  inputName = "contigs",
  outputName = "abricate",
  db = "ncbi",
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

  Name of ABRicate database to use (string)

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
