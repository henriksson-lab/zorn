# Get a Bascet binary from the target from a locally built Bascet repository

Get a Bascet binary from the target from a locally built Bascet
repository

## Usage

``` r
getBascetDevDir(
  devdir,
  tempdir = NULL,
  logLevel = "info",
  targetType = "release"
)
```

## Arguments

- tempdir:

  Default is to create a directory for temporary files in the current
  directory. Place it on a fast disk if possible

- logLevel:

  Log level for the Bascet instance (e.g. "info", "debug", "warn")

- targetType:

  What target type to load

- storeAt:

  Directory to store the container in. Default is current directory but
  it is likely better to provide a single systems level directory
  \#################### TODO

## Value

A Bascet instance
