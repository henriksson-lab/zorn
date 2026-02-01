# Filter an alignment (BAM-file).

This is typically used to either remove host DNA, or keep reads mapping
to a known reference. If the BAM-file has paired reads then BOTH reads
need to be mapped (flag 0x2); otherwise (flag 0x4)

## Usage

``` r
BascetFilterAlignment(
  bascetRoot,
  numThreads = NULL,
  inputName,
  outputName,
  keepMapped = FALSE,
  overwrite = FALSE,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- numThreads:

  Number of threads to use. Default is the maximum, taken from runner
  settings

- inputName:

  Name of input shards (BAM-file format)

- outputName:

  Name of output shards (BAM-file format)

- keepMapped:

  Keep the mapped reads (TRUE) or unmapped (FALSE)

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
