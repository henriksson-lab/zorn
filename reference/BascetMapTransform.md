# Transform data

This command enables

- subsetting to a list of cells

- converting between file formats

- merging shards

- dividing shards

## Usage

``` r
BascetMapTransform(
  bascetRoot,
  inputName,
  outputName,
  numDivide = 1,
  numMerge = 1,
  outFormat = "tirp.gz",
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

  Name of input shard

- outputName:

  Name of output shard

- numDivide:

  Must be \>=1; if more than 1, divide each input shard into number of
  outputs given here

- numMerge:

  Must be \>=1; if more than 1, merge this number of input shards into
  one

- outFormat:

  Extension for the output files

- includeCells:

  List of cells to include, or NULL if to include all

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
