# Aggregate data from previous Map call

Aggregate data from previous Map call

## Usage

``` r
BascetAggregateMap(
  bascetRoot,
  inputName,
  aggrFunction,
  includeCells = NULL,
  showProgress = TRUE,
  verbose = FALSE,
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- aggrFunction:

  Function to use for extracting simplified data for cell

- includeCells:

  Cells to aggregate

- showProgress:

  Show progress bar

- verbose:

  Show debug output

## Value

Aggregated data
