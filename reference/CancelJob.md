# Cancel a job

Cancel a job

## Usage

``` r
CancelJob(job)

# S4 method for class 'NoJob'
CancelJob(job)

# S4 method for class 'LocalJob'
CancelJob(job)

# S4 method for class 'SlurmJob'
CancelJob(job)
```

## Arguments

- job:

  A job object.

## Value

Backend-specific cancellation result, usually invisible.

## Methods (by class)

- `CancelJob(NoJob)`: No-op job method.

- `CancelJob(LocalJob)`: Local job method.

- `CancelJob(SlurmJob)`: SLURM job method.
