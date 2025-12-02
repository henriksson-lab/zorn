# From a Signac chromatin assay with fragments, for each cell, count how many reads per chromosome. This function directly returns an assay that can be added to a Seurat multimodal object

From a Signac chromatin assay with fragments, for each cell, count how
many reads per chromosome. This function directly returns an assay that
can be added to a Seurat multimodal object

## Usage

``` r
FragmentCountsPerChromAssay(bascetRoot, inputName = "fragments.1.tsv.gz")
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

## Value

A seurat object holding the counts
