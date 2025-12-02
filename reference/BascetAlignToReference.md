# Align from FASTQ, generate sorted and indexed BAM file

Align from FASTQ, generate sorted and indexed BAM file

## Usage

``` r
BascetAlignToReference(
  bascetRoot,
  useReference,
  numLocalThreads = NULL,
  inputName = "asfq",
  outputNameBAMunsorted = "unsorted_aligned",
  outputNameBAMsorted = "aligned",
  doSort = TRUE,
  overwrite = FALSE,
  aligner = c(NULL, "BWA", "STAR"),
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- useReference:

  Name of the BWA reference to use

- numLocalThreads:

  Number of threads to use for each runner. Default is the maximum,
  taken from runner settings

- inputName:

  Name of input shard

- outputNameBAMunsorted:

  Name of unsorted BAMs

- outputNameBAMsorted:

  Name of sorted BAMs (if generated)

- doSort:

  Whether to sort the output or not

- overwrite:

  Force overwriting of existing files. The default is to do nothing
  files exist

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance
