# Build count table from FASTQ reads and a list of selected kmers

Build count table from FASTQ reads and a list of selected kmers

## Usage

``` r
BascetQueryFq(
  bascetRoot,
  inputName = "filtered",
  outputName = "kmer_counts",
  useKMERs,
  maxReads = 1e+06,
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

- useKMERs:

  List of KMERs to query

- maxReads:

  The maximum number of reads per cell to sample

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
