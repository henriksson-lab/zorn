# Index a genome using BWA such that it can be used for alignment

FUTURE: could check if genome is indexed already

## Usage

``` r
BascetIndexGenomeBWA(
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
