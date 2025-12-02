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
