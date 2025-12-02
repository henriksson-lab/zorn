# Create a runner that submits jobs to SLURM

Create a runner that submits jobs to SLURM

## Usage

``` r
SlurmRunner(
  settings = NULL,
  ncpu = NULL,
  partition = NULL,
  account = NULL,
  time = NULL,
  prepend = NULL,
  mem = NULL,
  direct = NULL,
  deleteScript = NULL,
  benchmark = NULL,
  verbose = NULL,
  logTime = NULL
)
```

## Arguments

- settings:

  Default settings to override; can be NULL

- ncpu:

  Number of cores requested (SLURM -c)

- partition:

  Which partition to run the job on (SLURM -p)

- account:

  Which account to run the job on (SLURM -A)

- time:

  The time the job is allowed to run, e.g. "0-72:00:00" (SLURM -t)

- prepend:

  Something to prepend to the command. TODO seems not used. present in
  instance instead!!

- mem:

  Amount of main memory to reserve (SLURM –mem)

- direct:

  Run and get the result directly. FALSE implies asynchronous execution

- deleteScript:

  Delete job script after execution. Set to FALSE if you want to dissect
  it for debugging purposes

- benchmark:

  Enable logging of final CPU and memory usage

- verbose:

  Enable additional debug output

- logTime:

  Log execution time

## Value

A SLURM runner
