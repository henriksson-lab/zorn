# Aggregate data from Ariba This is a thin wrapper around BascetAggregateMap

Aggregate data from Ariba This is a thin wrapper around
BascetAggregateMap

## Usage

``` r
BascetAggregateAriba(
  bascetRoot,
  inputName = "ariba",
  verbose = NULL,
  includeCells = NULL,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- verbose:

  Show debug output

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

Aggregated data
