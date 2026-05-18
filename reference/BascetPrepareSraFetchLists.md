# Write sharded SRA accession lists for later import

`BascetPrepareSraFetchLists()` reads an SRA RunInfo table, orders runs
by cell name, and writes one `.sralist` file per zorn shard. These files
are intended as inputs to external download jobs that later feed
`bascet import-sra`.

## Usage

``` r
BascetPrepareSraFetchLists(
  runinfo,
  bascetRoot,
  outputName = "tofetch",
  numShards = NULL,
  runsPerShard = 1000,
  runColumn = "Run",
  cellColumn = NULL,
  overwrite = FALSE
)
```

## Arguments

- runinfo:

  SRA RunInfo data frame, or path to a RunInfo CSV file.

- bascetRoot:

  The root folder where all Bascets are stored.

- outputName:

  Name of output shard prefix.

- numShards:

  Number of `.sralist` files to create. If NULL, this is derived from
  `runsPerShard`.

- runsPerShard:

  Target number of SRA runs per `.sralist` when `numShards` is NULL.

- runColumn:

  Column containing SRA run accessions.

- cellColumn:

  Column used for alphabetic cell ordering.

- overwrite:

  Overwrite existing output files.

## Value

A data frame describing run-to-shard assignment.
