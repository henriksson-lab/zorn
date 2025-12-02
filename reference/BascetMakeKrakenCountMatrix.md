# Produce a count matrix of taxonomy IDs from KRAKEN output

Produce a count matrix of taxonomy IDs from KRAKEN output

## Usage

``` r
BascetMakeKrakenCountMatrix(
  bascetRoot,
  numLocalThreads = 1,
  inputName = "kraken_out",
  outputName = "kraken",
  overwrite = FALSE,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- numLocalThreads:

  Number of threads for KRAKEN to use

- inputName:

  Name of input shard (KRAKEN output)

- outputName:

  Name of output shard (h5 count matrix)

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
