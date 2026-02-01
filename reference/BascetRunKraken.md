# Run KRAKEN2 for each cell. Then produce a count matrix of taxonomy IDs from the output

Run KRAKEN2 for each cell. Then produce a count matrix of taxonomy IDs
from the output

## Usage

``` r
BascetRunKraken(
  bascetRoot,
  useKrakenDB = NULL,
  numThreads = NULL,
  inputName = "filtered",
  outputRawName = "kraken_raw",
  outputMatrixName = "kraken_mat",
  overwrite = FALSE,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- useKrakenDB:

  Path to KRAKEN2 database

- numThreads:

  Number of threads for one KRAKEN instance

- inputName:

  Name of input shard (FASTQ)

- outputRawName:

  Name of output shard (kraken raw output)

- outputMatrixName:

  Name of output shard (kraken count table data)

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
