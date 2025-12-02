# From aligned BAM file, compute counts per chromosome

From aligned BAM file, compute counts per chromosome

## Usage

``` r
BascetCountChrom(
  bascetRoot,
  inputName = "aligned",
  outputName = "chromcount",
  minMatching = 0,
  removeDuplicates = TRUE,
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

- minMatching:

  Disregard reads having fewer than specified matches, based on CIGAR
  string

- removeDuplicates:

  Deduplicate reads

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
