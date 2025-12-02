# From aligned BAM file, compute counts per feature

From aligned BAM file, compute counts per feature

## Usage

``` r
BascetCountFeature(
  bascetRoot,
  inputName = "aligned",
  outputName = "featurecount",
  gffFile,
  useFeature = "gene",
  attrGeneId = "gene_id",
  attrGeneName = "name",
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

- gffFile:

  GFF-like file to use for feature annotation

- useFeature:

  Feature type in the file to count

- attrGeneId:

  Attribute field to use this on for gene ID

- attrGeneName:

  Attribute field to use this on for gene name

- overwrite:

  Force overwriting of existing files. The default is to do nothing
  files exist

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A job, executing the counting
