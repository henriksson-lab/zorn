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

- overwrite:

  Force overwriting of existing files. The default is to do nothing
  files exist

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A job to be executed, or being executed, depending on runner settings
