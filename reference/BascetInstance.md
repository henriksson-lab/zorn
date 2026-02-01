# Create a new bascet instance. For advanced users only

Create a new bascet instance. For advanced users only

## Usage

``` r
BascetInstance(bin, tempdir, prependCmd = "", containerMem = "0B")
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

## Value

A Bascet instance
