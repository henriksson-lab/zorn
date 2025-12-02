# Detect metadata for raw input FASTQ files

\_R1 – common from illumina sequencer SRR\*\*\*\*\_1.fastq.gz – typical
from SRA

## Usage

``` r
DetectRawFileMeta(rawRoot, verbose = FALSE)
```

## Arguments

- rawRoot:

  Path to folder with FASTQ files

- verbose:

  Print additional information, primarily to help troubleshooting

## Value

A data frame with metadata for the raw input files

## Details

TODO Would be convenient to handle multiple samples, as sample1/xxx; in
this case, should prepend the sample name to the barcodes.

issue: when shardifying, good to keep info about what to merge. this
reduces the work plenty! could keep a list of which shards belong for
the next step
