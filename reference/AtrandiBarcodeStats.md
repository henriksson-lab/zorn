# Given a Bascet, produce a matrix showing for each combinatorial barcode, how many times it occurs across the cells. Presented as a 96-well plate matrix

This command assumes that cells are named as follows:
well1_well2_well3_well4, where e.g. well1 is in format G12

## Usage

``` r
AtrandiBarcodeStats(
  bascetRoot,
  inputName = "debarcoded",
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard (should be debarcoded reads)

- bascetInstance:

  A Bascet instance

## Value

Matrix showing coverage of each barcode
