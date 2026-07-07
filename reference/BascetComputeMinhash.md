# Compute minhashes for each cell.

Runs the native, streaming `bascet minhash-fq` command directly on each
TIRP shard (one array job task per shard). Unlike the old
`BascetMapCell` based path, this never materialises a whole cell's reads
in memory, so it is safe for wildly uneven (e.g. MDA) cells.

## Usage

``` r
BascetComputeMinhash(
  bascetRoot,
  inputName = "filtered",
  outputName = "minhash",
  kmerSize = 31,
  numMinhash = 1000,
  numThreads = NULL,
  overwrite = FALSE,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard (a TIRP file)

- outputName:

  Name of output shard

- kmerSize:

  The KMER size for the hashing (1..32)

- numMinhash:

  Number of minhashes (features) kept per cell

- numThreads:

  Number of minhash worker threads. Default is the runner CPU count

- overwrite:

  Force overwriting of existing files. The default is to do nothing if
  files exist

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A job to be executed, or being executed, depending on runner settings
