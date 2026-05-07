# Figure out if a BAM-file is a paired alignment or not

Figure out if a BAM-file is a paired alignment or not

## Usage

``` r
isBamPairedAlignment(
  fname,
  bascetInstance = GetDefaultBascetInstance(),
  maxRecords = 10000
)
```

## Arguments

- fname:

  Name of BAM-file

- bascetInstance:

  Deprecated; ignored

- maxRecords:

  Maximum number of BAM records to scan

## Value

TRUE if the BAM-file is a paired alignment
