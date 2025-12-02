# Callback function for aggregating min-hashes for each cell. To be called from BascetAggregateMap

Callback function for aggregating min-hashes for each cell. To be called
from BascetAggregateMap

## Usage

``` r
aggr.minhash(bascetFile, cellID, bascetInstance)
```

## Arguments

- bascetFile:

  Bascet file handle

- cellID:

  ID of cell to process

- bascetInstance:

  A Bascet instance

## Value

Minhash data (minhash.txt) for each cell
