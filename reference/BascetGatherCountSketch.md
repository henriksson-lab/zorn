# Gather all count sketches into a single count sketch matrix

Gather all count sketches into a single count sketch matrix

## Usage

``` r
BascetGatherCountSketch(
  bascetRoot,
  inputName = "countsketch",
  outputName = "countsketch_mat.csv",
  includeCells = NULL,
  kmerSize = 31,
  sketchSize = 4096,
  overwrite = FALSE,
  numLocalThreads = NULL,
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

  Name of output file

- includeCells:

  List of cells to process

- kmerSize:

  KMER length to use

- sketchSize:

  Size of the count sketch. Must be a power of two

- overwrite:

  Force overwriting of existing files. The default is to do nothing
  files exist

- numLocalThreads:

  Number of threads to use per job. Default is the number from the
  runner

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

  TODO produce a binary file format instead; gather files upon loading?

## Value

A job to be executed, or being executed, depending on runner settings
