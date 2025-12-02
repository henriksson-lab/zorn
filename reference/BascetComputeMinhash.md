# Compute minhashes for each cell. This is a thin wrapper around BascetMapCell

Compute minhashes for each cell. This is a thin wrapper around
BascetMapCell

## Usage

``` r
BascetComputeMinhash(
  bascetRoot,
  inputName = "filtered",
  outputName = "minhash",
  overwrite = FALSE,
  maxReads = 1e+05,
  kmerSize = 31,
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

- overwrite:

  Force overwriting of existing files. The default is to do nothing
  files exist

- maxReads:

  The maximum number of reads per cell to sample

- kmerSize:

  The KMER size for the hashing

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A job to be executed, or being executed, depending on runner settings
