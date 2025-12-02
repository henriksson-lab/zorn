# De novo assembly

First set up your Zorn/Bascet work directory as before.

Assembly can then performed using SKESA, via the MapCell system:

[(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r
### Assemble all genomes
BascetMapCell(
  bascetRoot,
  withfunction = "_skesa",
  inputName = "filtered",
  outputName = "skesa"
)
```

This produce contigs for each cell, all grouped together in “skesa”.
Note that both skesa and spades only assembles cells for which there are
sufficient number of reads.

``` r
#TODO extract an individual genome
```

``` r
#TODO QUAST for QC etc
```
