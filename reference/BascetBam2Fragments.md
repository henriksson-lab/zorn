# Take aligned BAM file and produce Fragments.tsv.gz, compatible with Signac ATAC-seq style analysis

Take aligned BAM file and produce Fragments.tsv.gz, compatible with
Signac ATAC-seq style analysis

## Usage

``` r
BascetBam2Fragments(
  bascetRoot,
  inputName = "aligned",
  outputName = "fragments",
  overwrite = FALSE,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- outputName:

  Name of output shard

- overwrite:

  Force overwriting of existing files. The default is to do nothing
  files exist

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A job, producing a type of Fragments.tsv.gz
