# Prepare to shard reads by collecting statistics about each barcode, and filtering out cells with few reads

Prepare to shard reads by collecting statistics about each barcode, and
filtering out cells with few reads

## Usage

``` r
PrepareSharding(
  bascetRoot,
  inputName = "debarcoded",
  minQuantile = 0.5,
  bascetInstance = GetDefaultBascetInstance(),
  verbose = TRUE
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- minQuantile:

  Read count-based cutoff for inclusion in final shards

- bascetInstance:

  A Bascet instance

- verbose:

  Print additional information, primarily to help troubleshooting

## Value

Statistics about the debarcoded reads
