# Index a genome using STAR such that it can be used for alignment

Index a genome using STAR such that it can be used for alignment

## Usage

``` r
BascetIndexGenomeSTAR(
  fastaFile,
  gtfFile,
  outDir,
  numLocalThreads = NULL,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- fastaFile:

  Name of FASTA file holding genome sequence

- gtfFile:

  GFF file holding genome annotation

- outDir:

  A directory in which to store the index. This directory will be
  created

- numLocalThreads:

  The number of threads to use for the STAR index, for each runner.
  Default is the maximum, taken from runner settings

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A job to be executed, or being executed, depending on runner settings
