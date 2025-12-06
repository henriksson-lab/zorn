# This creates a job object, linking to a running command. Mainly used for development but can be used in case a Zorn session died and you want to create a new monitor

This creates a job object, linking to a running command. Mainly used for
development but can be used in case a Zorn session died and you want to
create a new monitor

## Usage

``` r
createSlurmJobFromExisting(pid, arraysize)
```

## Arguments

- pid:

  Process ID

- arraysize:

  Size of array (must be known)

## Value

A SLURM job object
