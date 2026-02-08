# Convert data to Bascet-FASTQ

Convert data to Bascet-FASTQ

## Usage

``` r
BascetToFastq(
  bascetRoot,
  inputName,
  outputName,
  numLocalThreads = NULL,
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

  Number of threads to use per job. Default is the number from the
  runner

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
