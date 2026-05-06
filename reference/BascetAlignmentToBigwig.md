# Generate a bigwig out of all reads in a sorted BAM. Note that the caller is responsible for sorting the BAM first

This function is mainly for QC purposes.

## Usage

``` r
BascetAlignmentToBigwig(
  bascetRoot,
  inputName = "aligned_pos",
  outputName = "pileup",
  overwrite = FALSE,
  numThreads = NULL,
  totalMem = NULL,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard (BAM-file)

- outputName:

  Name of output shard (BIGWIG-file)

- overwrite:

  Force overwriting of existing files. The default is to do nothing
  files exist

- numThreads:

  Number of threads for BAM decompression and BigWig writing

- totalMem:

  Total memory to allocate

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A runner job (details depends on runner)
