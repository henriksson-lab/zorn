# Detect metadata for raw input FASTQ files

\_R1 – common from illumina sequencer SRR\*\*\*\*\_1.fastq.gz – typical
from SRA

## Usage

``` r
bascetCheckOverwriteOutput(outputFiles, overwrite)
```

## Arguments

- outputFiles:

  Files that we expect to exist

- overwrite:

  Files that we expect to exist

## Value

TRUE if ok to proceed

## Details

TODO Would be convenient to handle multiple samples, as sample1/xxx; in
this case, should prepend the sample name to the barcodes.

issue: when shardifying, good to keep info about what to merge. this
reduces the work plenty! could keep a list of which shards belong for
the next step
