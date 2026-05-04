# Run FASTQC on reads of all cells.

This uses the bascet integrated fastqc command rather than the old
mapcell system.

## Usage

``` r
BascetMapCellFASTQC(
  bascetRoot,
  inputName = "filtered",
  outputName = "fastqc",
  numThreads = NULL,
  numThreadsRead = NULL,
  nogroup = FALSE,
  expgroup = FALSE,
  kmerSize = 7,
  minLength = 0,
  dupLength = 50,
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

- numThreads:

  Total thread budget. Defaults to the runner CPU count

- numThreadsRead:

  Threads used by the TIRP reader. If NULL, use the CLI default

- nogroup:

  Do not group bases in the FastQC per-base modules

- expgroup:

  Use exponential base grouping in the FastQC per-base modules

- kmerSize:

  K-mer size for FastQC k-mer content

- minLength:

  Minimum sequence length to include

- dupLength:

  Length to truncate sequences for duplication detection

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
