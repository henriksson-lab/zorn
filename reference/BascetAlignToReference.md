# Align from FASTQ, generate sorted and indexed BAM file

Align from FASTQ, generate sorted and indexed BAM file

## Usage

``` r
BascetAlignToReference(
  bascetRoot,
  useReference,
  numThreads = NULL,
  totalMem = NULL,
  inputName = "filtered",
  outputName = "aligned",
  outputNameBAMcell = NULL,
  outputNameBAMpos = NULL,
  overwrite = FALSE,
  aligner = c(NULL, "BWAMEM2", "STAR", "minimap2"),
  bwamem2BatchPairs = NULL,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- useReference:

  Name of the BWA reference to use

- numThreads:

  Number of threads to use for each runner. Default is the maximum,
  taken from runner settings

- totalMem:

  Total memory to allocate

- inputName:

  Name of input shard

- outputName:

  Output shard name prefix. Defaults to `<outputName>_cell` for the
  unsorted/cell BAM and `<outputName>_pos` for the position-sorted BAM.

- outputNameBAMcell:

  Name of cell-sorted BAMs. If NULL, derived from `outputName`.

- outputNameBAMpos:

  Name of pos-sorted BAMs. If NULL, derived from `outputName`.

- overwrite:

  Force overwriting of existing files. The default is to do nothing
  files exist

- aligner:

  Which aligner to use: "BWAMEM2", "STAR", or "minimap2"

- bwamem2BatchPairs:

  For aligner "BWAMEM2" only: max read pairs per alignment batch. Bounds
  per-batch peak memory (bwa-mem2 materialises a whole batch's alignment
  scratch + SAM before returning). If NULL (default), bascet uses its
  built-in default (100000). Lower it on memory-constrained nodes, raise
  it for throughput.

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance
