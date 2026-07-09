# Wait for a job to finish

Wait for a job to finish

## Usage

``` r
WaitForJob(job)

# S4 method for class 'NoJob'
WaitForJob(job)

# S4 method for class 'SlurmJob'
WaitForJob(job)
```

## Arguments

- job:

  A job object.

## Value

Invisibly returns when the job has completed.

## Methods (by class)

- `WaitForJob(NoJob)`: No-op job method.

- `WaitForJob(SlurmJob)`: SLURM job method.
