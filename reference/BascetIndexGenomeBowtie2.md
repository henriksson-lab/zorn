# Index a genome using Bowtie2 such that it can be used for alignment

Index a genome using Bowtie2 such that it can be used for alignment

## Usage

``` r
BascetIndexGenomeBowtie2(
  genomeFile,
  numThreads = NULL,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- genomeFile:

  Name of FASTA file holding genome sequence

- numThreads:

  Number of threads to use. Default is the maximum, taken from runner
  settings

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance
