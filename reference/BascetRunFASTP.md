# Run FASTP for each cell. Input must be in FASTQ file format

Run FASTP for each cell. Input must be in FASTQ file format

## Usage

``` r
BascetRunFASTP(
  bascetRoot,
  inputName = "asfq",
  outputName = "fastp",
  numLocalThreads,
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

- numLocalThreads:

  Number of threads to use for FASTP. TODO: autodetect? take from
  runner? (fastp has no good defaults)

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
