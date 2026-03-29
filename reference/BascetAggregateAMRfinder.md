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

- includeCells:

  Character vector of cell names to include, or NULL for all cells

- getColumn:

  Column name from AMRfinder output to use as feature (string)

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

Aggregated data
