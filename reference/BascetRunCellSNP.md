# Align from FASTQ, generate sorted and indexed BAM file

TODO Should be managed by Bascet mapshard system, with automatic input
conversion. unaligned file should be made temp and removed

## Usage

``` r
BascetRunCellSNP(
  bascetRoot,
  inputName = "aligned",
  outputName = "cellsnp",
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

  Number of threads to use for each runner. Default is the maximum,
  taken from runner settings

- overwrite:

  Force overwriting of existing files. The default is to do nothing
  files exist

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance
