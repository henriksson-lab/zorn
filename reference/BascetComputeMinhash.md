# Compute minhashes for each cell. This is a thin wrapper around BascetMapCell

Compute minhashes for each cell. This is a thin wrapper around
BascetMapCell

## Usage

``` r
BascetComputeMinhash(
  bascetRoot,
  inputName = "filtered",
  outputName = "minhash",
  maxReads = 1e+05,
  kmerSize = 31,
  ...
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- outputName:

  Name of output shard

- maxReads:

  The maximum number of reads per cell to sample

- kmerSize:

  The KMER size for the hashing

- ...:

  Additional arguments passed to
  [`BascetMapCell`](https://henriksson-lab.github.io/zorn/reference/BascetMapCell.md)

## Value

A job to be executed, or being executed, depending on runner settings

## See also

[`BascetMapCell`](https://henriksson-lab.github.io/zorn/reference/BascetMapCell.md)
