# Downloading NCBI reference genomes

## When to use this

Zorn/Bascet can download public NCBI assemblies and pack them directly
into Bascet-ZIP shards. Each assembly is treated as one cell, and each
cell contains one file

## Download the metadata directory

This downloads the RefSeq bacterial assembly summary from NCBI and
caches it under the user’s cache directory:

``` r

metadataFile <- BascetDownloadNcbiGenomeMetadata(
  db = "refseq",
  group = "bacteria"
)

assemblies <- BascetReadNcbiGenomeMetadata(metadataFile)
dim(assemblies)
```

## Pick genomes

The table can be treated like a data.frame in R, and filtered to get the
genomes you are interested in. This is just a random example of
filtering:

``` r

selected <- assemblies[
  assemblies$version_status == "latest" &
    (is.na(assemblies$excluded_from_refseq) | assemblies$excluded_from_refseq == "") &
    assemblies$refseq_category %in% c("reference genome", "representative genome") &
    !is.na(assemblies$ftp_path) &
    nzchar(assemblies$ftp_path),
]
```

## Download NCBI genomes

We will run several downloads in parallel, each creating one shard. To
do this, we first pepare input lists for each shard:

``` r

BascetPrepareNcbiGenomeFetchLists(
  bascetRoot = bascetRoot,
  selected
)
```

Then start the download jobs:

``` r

BascetDownloadNcbiGenomes(
  bascetRoot = bascetRoot
)
```

By default, Bascet writes ZIP shards. To write TIRP shards instead, set
the suggested output extension:

``` r

BascetDownloadNcbiGenomes(
  bascetRoot = bascetRoot,
  inputName = "tofetch",
  outputName = "filtered",
  outFormat = "tirp.gz"
)
```

## NCBI load vs SLURM

The Bascet downloader uses a global download-start rate limit inside
each job. The default in the Zorn wrapper is conservative to avoid
problems. Note that if you launch many array jobs at once (SLURM), the
total rate is multiplied by the number of simultaneously running jobs.
If this is a problem depends on things like the IP of the compute nodes
etc.
