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
  outputNameBAMcell = "aligned_cell",
  outputNameBAMpos = "aligned_pos",
  overwrite = FALSE,
  aligner = c(NULL, "BWA", "bowtie2", "STAR"),
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

- outputNameBAMcell:

  Name of cell-sorted BAMs

- outputNameBAMpos:

  Name of pos-sorted BAMs (if generated)

- overwrite:

  Force overwriting of existing files. The default is to do nothing
  files exist

- aligner:

  Which aligner to use: "BWA", "bowtie2", or "STAR"

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance
