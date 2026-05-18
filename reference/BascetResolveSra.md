# Resolve SRA inputs to a canonical run table

Accepts direct SRR accessions, an `.sralist`, an existing RunInfo
CSV/data frame, or an NCBI SRA query. The returned table always contains
`Run` and `cell`, with `cell` defaulting to `Run` to avoid collisions.
Other RunInfo columns are retained as a mapping table for later
renaming.

## Usage

``` r
BascetResolveSra(
  runs = NULL,
  sralist = NULL,
  runinfo = NULL,
  query = NULL,
  bioproject = NULL,
  terms = NULL,
  pageSize = 500,
  sleep = 1.2
)
```

## Arguments

- runs:

  Character vector of SRA run accessions.

- sralist:

  Path to a file with one SRA run accession per line.

- runinfo:

  SRA RunInfo data frame, or path to a RunInfo CSV file.

- query:

  NCBI SRA query passed to `rentrez`.

- bioproject:

  BioProject accession. Used to build a query when `query` is not
  provided.

- terms:

  Additional NCBI SRA query terms joined with `AND`.

- pageSize:

  Number of records fetched per E-utilities page.

- sleep:

  Seconds to sleep between paged E-utilities fetches.

## Value

A canonical SRA run table.
