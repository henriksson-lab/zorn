# Genome annotation tools

## QUAST

First run QUAST on all the cells: [(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r
BascetMapCell(
  bascetRoot,
  withfunction = "_quast",
  inputName = "skesa",  #or other source of contigs
  outputName = "quast"
)
```

Then aggregate the results for visualization. This example caches the
result to speed up reloading; this is optional

``` r
quast_aggr <- BascetCacheComputation(bascetRoot,"cache_quast",MapListAsDataFrame(BascetAggregateMap(
  bascetRoot,
  "quast",
  aggr.quast
)))
```

## Abricate

First run Abricate on all the cells. The NCBI database is used by
default. See ListDatabaseAbricate() for a list of other databases

[(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r
BascetMapCellAbricate(
  bascetRoot,
  inputName = "skesa"  #or other source of contigs
)
```

Then aggregate the results for visualization. This example caches the
result to speed up reloading; this is optional

``` r
quast_aggr <- BascetCacheComputation(bascetRoot,"cache_abricate",MapListAsDataFrame(BascetAggregateMap(
  bascetRoot,
  "abricate",
  aggr.abricate
)))
```

## Bakta

First download a database:

``` r
DownloadDatabaseBakta(
  dbdir="~/bakta",  #create directory before running command
  dbtype="light"
)
```

You can run then run Bakta on all cells: [(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r
### ....
BascetMapCellBakta(
  bascetRoot,
  db = "~/bakta",
  inputName = "skesa",  #or other source of contigs
)
```

Then aggregate the results for visualization. This example caches the
result to speed up reloading; this is optional

``` r
quast_aggr <- BascetCacheComputation(bascetRoot,"cache_bakta",MapListAsDataFrame(BascetAggregateMap(
  bascetRoot,
  "bakta",
  aggr.bakta
)))
```
