# Generate a bigwig out of all reads in a sorted BAM. Note that the caller is responsible for sorting the BAM first

This function is mainly for QC purposes. It uses bamCoverage from
deepTools apt install python3-deeptools

## Usage

``` r
BascetAlignmentToBigwig(
  bascetRoot,
  inputName = "aligned",
  outputName = "pileup",
  overwrite = FALSE,
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

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A runner job (details depends on runner)
