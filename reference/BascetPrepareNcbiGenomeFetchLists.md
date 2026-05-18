# Write sharded NCBI genome download inputs

`BascetPrepareNcbiGenomeFetchLists()` filters/normalizes an NCBI
assembly metadata table and writes one TSV input shard per Bascet
download job. The generated shards contain `cell_id`,
`assembly_accession`, and `ftp_path`.

## Usage

``` r
BascetPrepareNcbiGenomeFetchLists(
  assemblies,
  bascetRoot,
  outputName = "tofetch",
  numShards = NULL,
  genomesPerShard = 1000,
  latest = TRUE,
  excludeFromRefseq = TRUE,
  assemblyLevel = NULL,
  refseqCategory = NULL,
  overwrite = FALSE
)
```

## Arguments

- assemblies:

  Assembly metadata data frame, or path to `assembly_summary.txt`.

- bascetRoot:

  The root folder where all Bascets are stored.

- outputName:

  Name of output shard prefix.

- numShards:

  Number of TSV shards to create. If NULL, derived from
  `genomesPerShard`.

- genomesPerShard:

  Target number of genomes per shard when `numShards` is NULL.

- latest:

  Keep only rows where `version_status == "latest"`.

- excludeFromRefseq:

  Keep only rows where `excluded_from_refseq` is empty.

- assemblyLevel:

  Optional assembly levels to keep.

- refseqCategory:

  Optional RefSeq categories to keep.

- overwrite:

  Overwrite existing output files.

## Value

A data frame describing genome-to-shard assignment.
