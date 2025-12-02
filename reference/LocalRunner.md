# Create new local runner instance

Create new local runner instance

## Usage

``` r
LocalRunner(ncpu = NULL, direct = TRUE, showScript = FALSE)
```

## Arguments

- ncpu:

  Number of CPU cores to use. By default, will attempt to detect and use
  all of them

- direct:

  Run jobs synchronously

- showScript:

  Show the script to run, for debugging purposes

## Value

A runner instance
