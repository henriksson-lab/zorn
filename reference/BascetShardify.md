# Take debarcoded reads, merge them, and split them into suitable numbers of shards.

The reads from one cell is guaranteed to only be present in a single
shard. This makes parallel processing simple as each shard can be
processed on a separate computer. Using more shards means that more
computers can process the data in parallel. However, if you perform all
the calculations on a single computer, having more than one shard will
not result in a speedup. This option is only relevant when using a
cluster of compute nodes.

## Usage

``` r
BascetShardify(
  debstat,
  numOutputShards = 1,
  outputName = "filtered",
  overwrite = FALSE,
  bufferSize = 1500,
  pageSize = 32,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- debstat:

  Plan for sharding provided by PrepareSharding

- numOutputShards:

  How many shards to generate /for each input prefix/

- outputName:

  Name of the output file: Properly sharded debarcoded reads

- overwrite:

  Force overwriting of existing files. The default is to do nothing
  files exist

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A runner job (details depends on runner)
