# Run GECCO on contigs of all cells.

This uses the bascet integrated gecco command rather than the old
mapcell system.

## Usage

``` r
BascetMapCellGECCO(
  bascetRoot,
  inputName = "contigs",
  outputName = "gecco",
  numThreads = NULL,
  dataDir = NULL,
  threshold = 0.8,
  cds = 3,
  noMask = FALSE,
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

- numThreads:

  Total thread budget. Defaults to the runner CPU count

- dataDir:

  GECCO data directory containing HMM, CRF model, and InterPro files

- threshold:

  Minimum probability for cluster membership

- cds:

  Minimum number of annotated CDS in a cluster

- noMask:

  Do not mask ambiguous nucleotides during gene prediction

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
