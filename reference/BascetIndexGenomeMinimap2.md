# Index a genome using minimap2 such that it can be used for alignment

Index a genome using minimap2 such that it can be used for alignment

## Usage

``` r
BascetIndexGenomeMinimap2(
  fastaFile,
  indexFile,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- fastaFile:

  Name of FASTA file holding genome sequence

- indexFile:

  Path to write the resulting .mmi index

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A job to be executed, or being executed, depending on runner settings
