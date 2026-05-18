# Download sharded NCBI genome inputs into Bascet files

Runs `bascet ncbi-genome-download` once per `.ncbilist` shard produced
by
[`BascetPrepareNcbiGenomeFetchLists()`](https://henriksson-lab.github.io/zorn/reference/BascetPrepareNcbiGenomeFetchLists.md).

## Usage

``` r
BascetDownloadNcbiGenomes(
  bascetRoot,
  inputName = "tofetch",
  outputName = "contigs",
  outFormat = "zip",
  threads = NULL,
  downloadStartsPerSecond = 2,
  queueSize = NULL,
  maxRetries = 5,
  keepTemp = FALSE,
  overwrite = FALSE,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored.

- inputName:

  Name of input `.ncbilist` shard prefix.

- outputName:

  Name of output shard prefix.

- outFormat:

  Output shard extension. Use `"zip"` to store `{assembly}/contigs.fa`
  entries, or `"tirp.gz"` to store contigs as R1 reads with dummy `F`
  quality scores and empty R2/UMI.

- threads:

  Parallel genome workers per Bascet job. Defaults to runner ncpu.

- downloadStartsPerSecond:

  Global per-job NCBI download start rate.

- queueSize:

  Optional completed-fragment queue size.

- maxRetries:

  Download retry count.

- keepTemp:

  Keep per-job temporary download and fragment files.

- overwrite:

  Whether to submit jobs when output files already exist.

- runner:

  Bascet runner.

- bascetInstance:

  Bascet instance.

## Value

A zorn job object.
