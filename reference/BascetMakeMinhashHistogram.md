# Gather all minhashes into a single histogram file

Gather all minhashes into a single histogram file

## Usage

``` r
BascetMakeMinhashHistogram(
  bascetRoot,
  inputName = "minhash",
  outputName = "minhash_hist.csv",
  includeCells = NULL,
  overwrite = FALSE,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard (should contain minhashes)

- outputName:

  Name of output file

- includeCells:

  List of cells to include, or NULL for all cells

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
