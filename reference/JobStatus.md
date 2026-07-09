# Get job status

Get job status

## Usage

``` r
JobStatus(job)

# S4 method for class 'NoJob'
JobStatus(job)

# S4 method for class 'LocalJob'
JobStatus(job)

# S4 method for class 'SlurmJob'
JobStatus(job)
```

## Arguments

- job:

  A job object.

## Value

A data frame or backend-specific status object.

## Methods (by class)

- `JobStatus(NoJob)`: No-op job method.

- `JobStatus(LocalJob)`: Local job method.

- `JobStatus(SlurmJob)`: SLURM job method.
