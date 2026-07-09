# Run a job

Submit a job script with a runner backend.

## Usage

``` r
RunJob(runner, jobname, bascetInstance, cmd, arraysize)

# S4 method for class 'NoRunner'
RunJob(runner, jobname, bascetInstance, cmd, arraysize)

# S4 method for class 'LocalRunner'
RunJob(runner, jobname, bascetInstance, cmd, arraysize)

# S4 method for class 'SlurmRunner'
RunJob(runner, jobname, bascetInstance, cmd, arraysize)
```

## Arguments

- runner:

  Runner object used to execute the job.

- jobname:

  Name for the submitted job.

- bascetInstance:

  Bascet instance that provides runtime settings.

- cmd:

  JobScript object to execute.

- arraysize:

  Number of array tasks.

## Value

A job object, or a no-op job for synchronous/no-runner backends.

## Methods (by class)

- `RunJob(NoRunner)`: No-op runner method.

- `RunJob(LocalRunner)`: Local runner method.

- `RunJob(SlurmRunner)`: SLURM runner method.
