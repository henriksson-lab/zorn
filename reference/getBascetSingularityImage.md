# Get a Bascet image (singularity or docker). It will be cached in the provided directory to avoid downloading it each the time the function is called

Get a Bascet image (singularity or docker). It will be cached in the
provided directory to avoid downloading it each the time the function is
called

## Usage

``` r
getBascetSingularityImage(storeAt = getwd(), tempdir = NULL, logLevel = "info")
```

## Arguments

- storeAt:

  Directory to store the container in. Default is current directory but
  it is likely better to provide a single systems level directory

- tempdir:

  Default is to create a directory for temporary files in the current
  directory. Place it on a fast disk if possible

## Value

A Bascet instance
