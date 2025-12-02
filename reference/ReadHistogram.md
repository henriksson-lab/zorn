# Read the count histogram associated with a Bascet. Not all Bascets have one, but it is typically produced after debarcoding

Read the count histogram associated with a Bascet. Not all Bascets have
one, but it is typically produced after debarcoding

## Usage

``` r
ReadHistogram(
  bascetRoot,
  inputName,
  bascetInstance = GetDefaultBascetInstance(),
  verbose = TRUE
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- bascetInstance:

  A Bascet instance

- verbose:

  Print additional information, primarily to help troubleshooting

## Value

Histogram as a data.frame
