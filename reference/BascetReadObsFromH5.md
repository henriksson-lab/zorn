# Read observation metadata from Bascet HDF5/AnnData files

Read observation metadata from Bascet HDF5/AnnData files

## Usage

``` r
BascetReadObsFromH5(bascetRoot, inputName, verbose = FALSE)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input HDF5 shard

- verbose:

  Print additional information, primarily to help troubleshooting

## Value

A data.frame with rownames set to cell names
