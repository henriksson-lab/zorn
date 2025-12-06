# Extract barcodes and trim input raw FASTQ

Extract barcodes and trim input raw FASTQ

## Usage

``` r
BascetGetRaw(
  bascetRoot,
  rawmeta,
  maxShardSize = "50g",
  outputName = "debarcoded",
  outputNameIncomplete = "incomplete_reads",
  chemistry = c("atrandi-wgs", "atrandi-rnaseq", "parse-bio"),
  subchemistry = NULL,
  barcodeTolerance = NULL,
  numLocalThreads = NULL,
  bufferSize = 4000,
  sortBufferSize = 4000,
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

- maxShardSize:

  Estimated maximum size of output shard. Can be set higher but as
  sorting is also performed during sharding, it can be overall more
  efficient to only do partial sorting during this command

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

  Number of threads to use per job. Default is the number from the
  runner

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A job to be executed, or being executed, depending on runner settings
