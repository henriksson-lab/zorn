# Create a new bascet instance. For advanced users only

Create a new bascet instance. For advanced users only

## Usage

``` r
BascetInstance(
  bin,
  tempdir,
  prependCmd = "",
  containerMem = "0B",
  logLevel = "info"
)
```

## Arguments

- bin:

  Name of the binary

- tempdir:

  Directory where to store temporary files

- prependCmd:

  Something to prepend to the command, to e.g. support container systems

- containerMem:

  Amount of memory used by the container itself

- logLevel:

  ...

## Value

A Bascet instance
