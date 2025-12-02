# Open a Bascet, prepare it for reading individual files

Note, the current code is based on pure R, but more efficient calls can
be made in the future. We thus advise against direct zip-file
manipulation and do not guarantee future support for this

## Usage

``` r
OpenBascet(bascetRoot, bascetName, bascetInstance = GetDefaultBascetInstance())
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- bascetName:

  Name of the bascet

- bascetInstance:

  A Bascet instance

## Value

A handle to a Bascet
