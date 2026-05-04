# Get a Bascet binary for the current platform It will be cached in the provided directory to avoid downloading it each the time the function is called

Get a Bascet binary for the current platform It will be cached in the
provided directory to avoid downloading it each the time the function is
called

## Usage

``` r
getBascetBinary(
  storeAt = NULL,
  tempdir = NULL,
  logLevel = "info",
  forceInstall = FALSE
)
```

## Arguments

- storeAt:

  Directory to store the binary in. If NULL, uses
  [`defaultBascetBinDir`](https://henriksson-lab.github.io/zorn/reference/defaultBascetBinDir.md)

- tempdir:

  Default is to create a directory for temporary files in the current
  directory. Place it on a fast disk if possible

- logLevel:

  Log level for the Bascet instance (e.g. "info", "debug", "warn")

- forceInstall:

  Force download of the Bascet binary even if a cached binary exists

## Value

A Bascet instance
