# Using SLURM

## Introduction

Zorn is designed to be a full replacement for the NextFlow job manager,
for the processing of single-cell genomics and microbial isolates. This
design is motivated by the Bascet file format: While NextFlow operates
on a file-by-file basis, Bascet keeps cells together in sharded objects.
These are natural minibatches of otherwise rather small tasks, reducing
the overhead of scheduling work for the (possibly many!) cells.

[SLURM](https://slurm.schedmd.com) is one of the most common job
managers out there, and Zorn thus works on top of SLURM by writing out
the appropriate Bascet commands. If you do not have SLURM on your
cluster, we are interested to hear what else people are using!

## Usage

To use SLURM, you have to get used to slightly different semantics than
if you run locally. First you need to set up a runner. The way we
recommend you to do it, these are only default values, such as which
partition to submit to, and the name of your account. You can also set a
permissible default timeout time:

``` r
## Set the new default runner. Note that the name should not be changed!
bascetRunner.default <- SlurmRunner(
  partition="shared", 
  account="name_of_your_project", 
  time="0-24:00:00"
)

## get a Bascet instance as usual; this is a separate concern from the runner
bascetInstance.default <- getBascetSingularityImage(storeAt="~/")
```

You then make calls to Bascet as in all other tutorials. The main
script, in charge of the control flow, should be run in
[tmux](https://github.com/tmux/tmux/wiki) or
[screen](https://linuxize.com/post/how-to-use-linux-screen/), on a login
node of your cluster.

You likely want to tune each SLURM job a bit in terms of how many cores
are used. This is the suggested method, where only relevant settings are
overridden on a per-job basis:

``` r
BascetGetRaw(
  bascetRoot,
  rawmeta,
  chemistry="atrandi-wgs",
  runner=SlurmRunner(
    bascetRunner.default,
    ncpu=5,
    mem="5g")
)
```

## SLURM settings

We are unable to give specific suggestions on cpu and memory as it
depends on how your SLURM instance is configured. As an example, Swedish
clusters are set up to give memory in proportion to the number of CPUs.

Bascet is written to avoid using large amounts of RAM when possible.
Exceptions are the use of downstream software such as KRAKEN (in
particular).

All operations are multithreaded and benefit from having more CPUs per
node. Speed is however limited by disk access speed. It it thus worth
investigating your log files after execution to see how much CPU/memory
you actually needed, to improve long-term running of Bascet.

## Advanced: Asynchronous execution

You can also set up SLURM jobs to run in the background, instead of
stalling your R session until done. This can primarily be useful if you
want to use the same R session for doing other work while waiting for
the job.

To set up this mode of running, you create the SLURM runner slightly
differently:

``` r
## Set the new default runner. Note that the name should not be changed!
bascetRunner.default <- SlurmRunner(
  ... #previous settings
  direct=FALSE ## NEW! this enables asynchronous mode
)

## get a Bascet instance as usual; this is a separate concern from the runner
bascetInstance.default <- getBascetSingularityImage(storeAt="~/")
```

Whenever you run a job, you now want to catch a handle to the job in a
variable:

``` r
my_job <- BascetGetRaw(  #note, store job in a variable
  bascetRoot,
  rawmeta,
  chemistry="atrandi-wgs",
  runner=SlurmRunner(
    bascetRunner.default,
    ncpu=5,
    mem="5g")
)
```

After running this command, you will immediately get the control back.
To see the status of a job, it is best to just run WaitForJob:

``` r
WaitForJob(my_job)
```

You can press ctrl+c to end the waiting loop early
