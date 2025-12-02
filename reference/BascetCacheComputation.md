# A wrapper to cache a computation. Put your function in as an argument, as R will only compute its value if needed. If the cache file exist, it will not be run again

A wrapper to cache a computation. Put your function in as an argument,
as R will only compute its value if needed. If the cache file exist, it
will not be run again

## Usage

``` r
BascetCacheComputation(bascetRoot, fname, value)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- fname:

  Name of the file to store the cache in. The extension .RDS is added
  automatically

- value:

  The value to be cached

## Value

The cached value
