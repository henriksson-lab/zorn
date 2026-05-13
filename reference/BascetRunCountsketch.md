# Gather all count sketches into a single count sketch matrix

Gather all count sketches into a single count sketch matrix

## Usage

``` r
BascetRunCountsketch(
  bascetRoot,
  inputName = "filtered",
  outputName = "countsketch_mat.feather",
  includeCells = NULL,
  kmerSize = 31,
  sketchSize = 4096,
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

- numThreads:

  Number of threads to use per job. Default is the number from the
  runner

- totalMem:

  Total memory to allocate

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A job to be executed, or being executed, depending on runner settings
