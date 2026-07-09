# Get job log

Get job log

## Usage

``` r
JobLog(job)

# S4 method for class 'LocalJob'
JobLog(job)

# S4 method for class 'SlurmJob'
JobLog(job)
```

## Arguments

- job:

  A job object.

## Value

Backend-specific log content, usually invisible if no log is available.

## Methods (by class)

- `JobLog(LocalJob)`: Local job method.

- `JobLog(SlurmJob)`: SLURM job method.
