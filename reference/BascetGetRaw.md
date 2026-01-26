# Extract barcodes and trim input raw FASTQ

Extract barcodes and trim input raw FASTQ

## Usage

``` r
BascetGetRaw(
  bascetRoot,
  rawmeta,
  maxShardSize = "50g",
  outputName = "debarcoded",
  chemistry = c("atrandi-wgs", "atrandi-rnaseq", "parse-bio"),
  barcodeTolerance = NULL,
  numThreads = NULL,
  numReadThreads = NULL,
  numDebarcodeThreads = NULL,
  numSortingThreads = NULL,
  numMergeSortingThreads = NULL,
  numWriteThreads = NULL,
  numCompressThreads = NULL,
  totalMem = NULL,
  streamBufferSize = NULL,
  sortBufferSize = NULL,
  compressBufferSize = NULL,
  compressRawBufferSize = NULL,
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

- chemistry:

  The type of data to be parsed

- barcodeTolerance:

  Optional: Number of mismatches allowed in the barcode for it to still
  be considered valid

- numThreads:

  Number of threads to use per job. Default is the number from the
  runner

- numReadThreads:

  Advanced setting: Number of threads for reading

- numDebarcodeThreads:

  Advanced setting: Number of threads for debarcoding

- numSortingThreads:

  Advanced setting: Number of threads for sorting, first phase

- numMergeSortingThreads:

  Advanced setting: Number of threads for sorting, second phase

- numWriteThreads:

  Advanced setting: Number of threads for writing

- numCompressThreads:

  Advanced setting: Number of threads for compressing

- totalMem:

  Total memory to allocate

- streamBufferSize:

  Advanced setting: Stream buffer size (fraction, given as e.g. "10%")

- sortBufferSize:

  Advanced setting: Sort buffer size (fraction, given as e.g. "10%")

- overwrite:

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

- pageBufferSize:

  Advanced setting: Page buffer size (fraction, given as e.g. "10%")

## Value

A job to be executed, or being executed, depending on runner settings
