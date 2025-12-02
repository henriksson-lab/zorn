# Take debarcoded reads and split them into suitable numbers of shards.

The reads from one cell is guaranteed to only be present in a single
shard. This makes parallel processing simple as each shard can be
processed on a separate computer. Using more shards means that more
computers can be used.

## Usage

``` r
BascetShardifyOld(
  bascetRoot,
  inputName = "debarcoded",
  includeCells = NULL,
  numOutputShards = 1,
  outputName = "filtered",
  overwrite = FALSE,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- inputName:

  Name of input file: Debarcoded reads

- outputName:

  Name of the output file: Properly sharded debarcoded reads

## Value

A job to be executed, or being executed, depending on runner settings

## Details

If you perform all the calculations on a single computer, having more
than one shard will not result in a speedup. This option is only
relevant when using a cluster of compute nodes.

TODO if we have multiple input samples, is there a way to group them?
otherwise we will be reading more input files than needed. that said, if
we got an index, so if list of cells specified, it is possible to
quickly figure out out if a file is needed at all for a merge

TODO Figuring out if a file is needed can be done at "planning" (Zorn)
stage

TODO seems faster to have a single merger that writes multiple output
files if cell list is not provided. if the overhead is accepted then
read all input files and discard cells on the fly
