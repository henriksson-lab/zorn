# Index a genome using BWA-MEM2 such that it can be used for alignment

TODO: could check if genome is indexed already

## Usage

``` r
BascetIndexGenomeBWAMEM2(
  genomeFile,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- genomeFile:

  Name of FASTA file holding genome sequence

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance
