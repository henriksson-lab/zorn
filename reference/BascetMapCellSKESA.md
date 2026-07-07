# Run integrated SKESA on reads from each cell.

Run integrated SKESA on reads from each cell.

## Usage

``` r
BascetMapCellSKESA(
  bascetRoot,
  inputName = "filtered",
  outputName = "contigs",
  numThreads = NULL,
  numSkesaWorkers = NULL,
  numSkesaCores = NULL,
  numThreadsRead = NULL,
  totalMem = NULL,
  kmer = 21,
  maxKmer = 0,
  steps = 11,
  minCount = NULL,
  maxKmerCount = 10,
  vectorPercent = 0.05,
  insertSize = 0,
  fraction = 0.1,
  maxSnpLen = 150,
  minContig = 50,
  allowSnps = FALSE,
  forceSingleEnds = FALSE,
  singlePassCounter = FALSE,
  maxReadsPerCell = 0,
  overwrite = FALSE,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- outputName:

  Name of output shard

- numThreads:

  Total thread budget. Defaults to the runner CPU count

- numSkesaWorkers:

  (Advanced) Number of cells to assemble concurrently

- numSkesaCores:

  (Advanced) Number of cores to give each SKESA assembly

- numThreadsRead:

  (Advanced) Threads used by the TIRP reader. If NULL, use the CLI
  default

- totalMem:

  Total memory to allocate

- kmer:

  Minimal k-mer length for assembly

- maxKmer:

  Maximal k-mer length for assembly. 0 means auto

- steps:

  Number of assembly iterations from minimal to maximal k-mer length

- minCount:

  Minimal count for k-mers retained. NULL (default) lets skesa
  auto-estimate it from coverage, raising it for high-coverage cells

- maxKmerCount:

  Maximum k-mer count for fork tie-breaking

- vectorPercent:

  Percentage of reads containing 19-mer for adapter detection. 1.0
  disables

- insertSize:

  Expected insert size for paired reads. 0 means auto

- fraction:

  Maximum noise to signal ratio acceptable for extension

- maxSnpLen:

  Maximal SNP length

- minContig:

  Minimal contig length reported in output

- allowSnps:

  Allow additional step for SNP discovery

- forceSingleEnds:

  Do not use paired-end information

- singlePassCounter:

  Use the legacy single-pass k-mer counter. Faster for small cells but
  can exceed the memory budget for high-coverage cells, so it is off by
  default

- maxReadsPerCell:

  Maximum read pairs per cell fed to assembly. 0 (default) disables the
  cap. When a cell exceeds this, only the first N read pairs encountered
  in the file are used (no random subsampling). Bounds memory for
  pathological high-read cells

- overwrite:

  Force overwriting of existing files. The default is to do nothing
  files exist

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A job to be executed, or being executed, depending on runner settings
