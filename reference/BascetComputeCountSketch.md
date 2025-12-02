# Compute count sketch for each cell. This is a thin wrapper around BascetMapCell

Compute count sketch for each cell. This is a thin wrapper around
BascetMapCell

## Usage

``` r
BascetComputeCountSketch(
  bascetRoot,
  inputName = "filtered",
  outputName = "countsketch",
  overwrite = FALSE,
  maxReads = 1e+05,
  kmerSize = 31,
  sketchSize = 5000,
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

  Max number of reads to sample per cell

- kmerSize:

  Size of the KMER to hash

- sketchSize:

  Number of dimensions to reduce to

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A job to be executed, or being executed, depending on runner settings
