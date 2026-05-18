# Download sharded SRA run lists into Bascet TIRP files

Runs `bascet import-sra` once per `.sralist` shard. The Bascet command
calls SRA Toolkit and writes indexed TIRP output directly.

## Usage

``` r
BascetDownloadSraRuns(
  bascetRoot,
  inputName = "tofetch",
  outputName = "filtered",
  runinfo = NULL,
  threads = NULL,
  runsAhead = 10,
  sraWorkers = NULL,
  overwrite = FALSE,
  keepTemp = FALSE,
  prefetch = "prefetch",
  fasterqDump = "fasterq-dump",
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored.

- inputName:

  Name of input `.sralist` shard prefix.

- outputName:

  Name of output TIRP shard prefix.

- runinfo:

  Optional RunInfo CSV. Defaults to `<inputName>.runinfo.csv` if present
  in `bascetRoot`.

- threads:

  Threads passed to `fasterq-dump`. Defaults to runner ncpu.

- runsAhead:

  Maximum number of completed `fasterq-dump` outputs buffered ahead of
  TIRP writing.

- sraWorkers:

  Number of concurrent SRA Toolkit workers. If NULL, Bascet uses
  `threads`.

- overwrite:

  Whether to submit jobs when output files already exist.

- keepTemp:

  Keep per-run SRA/FASTQ scratch files.

- prefetch:

  Path to SRA Toolkit `prefetch` executable.

- fasterqDump:

  Path to SRA Toolkit `fasterq-dump` executable.

- runner:

  Bascet runner.

- bascetInstance:

  Bascet instance.

## Value

A zorn job object.
