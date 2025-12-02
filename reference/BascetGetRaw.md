# Generate BAM with barcodes from input raw FASTQ

Generate BAM with barcodes from input raw FASTQ

## Usage

``` r
BascetGetRaw(
  bascetRoot,
  rawmeta,
  outputName = "debarcoded",
  outputNameIncomplete = "incomplete_reads",
  chemistry = c("atrandi_wgs", "atrandi_rnaseq", "pb_rnaseq"),
  subchemistry = NULL,
  barcodeTolerance = NULL,
  numLocalThreads,
  overwrite = FALSE,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- rawmeta:

  Metadata for the raw FASTQ input files. See DetectRawFileMeta

- outputName:

  Name output files: Debarcoded reads

- outputNameIncomplete:

  Name of output files: Reads that could not be parsed

- chemistry:

  The type of data to be parsed

- barcodeTolerance:

  Optional: Number of mismatches allowed in the barcode for it to still
  be considered valid

- numLocalThreads:

  Number of threads to use per job

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A job to be executed, or being executed, depending on runner settings
