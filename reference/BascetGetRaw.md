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
  numReadThreads = NULL,
  numDebarcodeThreads = NULL,
  numSortingThreads = NULL,
  numWriteThreads = NULL,
  streamBufferSize = NULL,
  sortBufferSize = NULL,
  pageBufferSize = NULL,
  totalMem = NULL,
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

- numReadThreads:

  Number of threads for reading (advanced; parameter not checked)

- numDebarcodeThreads:

  Number of threads for debarcoding (advanced; parameter not checked)

- numSortingThreads:

  Number of threads for sorting (advanced; parameter not checked)

- numWriteThreads:

  Number of threads for writing (advanced; parameter not checked)

- streamBufferSize:

  Stream buffer size (advanced; parameter not checked)

- sortBufferSize:

  Sort buffer size (advanced; parameter not checked)

- pageBufferSize:

  Page buffer size (advanced; parameter not checked)

- totalMem:

  Total memory to allocate

- overwrite:

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A job to be executed, or being executed, depending on runner settings
