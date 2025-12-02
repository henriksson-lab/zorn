# Run KRAKEN2 for each cell

Run KRAKEN2 for each cell

## Usage

``` r
BascetRunKraken(
  bascetRoot,
  useKrakenDB = "/data/henlab/kraken/standard-8",
  numLocalThreads = NULL,
  inputName = "asfq",
  outputName = "kraken_out",
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

- numLocalThreads:

  Number of threads for one KRAKEN instance

- inputName:

  Name of input shard (FASTQ)

- outputName:

  Name of output shard (kraken data)

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
