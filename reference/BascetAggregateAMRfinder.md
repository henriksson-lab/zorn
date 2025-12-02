# Aggregate data from AMRfinder This is a thin wrapper around BascetAggregateMap

Aggregate data from AMRfinder This is a thin wrapper around
BascetAggregateMap

## Usage

``` r
BascetAggregateAMRfinder(
  bascetRoot,
  inputName = "AMRfinder",
  includeCells = NULL,
  getColumn = "Element.symbol",
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

- verbose:

  Show debug output

## Value

Aggregated data
