# Doublet removal

First set up Zorn/Bascet according to the [install
instructions](https://henriksson-lab.github.io/zorn/articles/install.md).
This tutorial assumes that you have already produced a Seurat object
using one of the count workflows, e.g. the [KRAKEN2
workflow](https://henriksson-lab.github.io/zorn/articles/kraken.md), the
[alignment-based
workflow](https://henriksson-lab.github.io/zorn/articles/alignment.md),
the [informative KMER
workflow](https://henriksson-lab.github.io/zorn/articles/kmer.md), or
the [count sketch
workflow](https://henriksson-lab.github.io/zorn/articles/countsketch.md).
You should have a count table of some sorts ready for analysis.

## When to remove doublets

A *doublet* in single-cell metagenomics is a barcode that captures DNA
from more than one cell. These can arise for at least two common
reasons:

1.  If more cells are collected than there are barcodes, you get
    **barcode collisions**. Droplet and split-pool methods rely on the
    cells being much fewer than the avaiable barcodes. By doing “Poisson
    loading” (diluting enough), doublets are typically kept \<2%; but
    they still occur and must be handled

2.  If cells stick together, the will also receive the same barcode.
    Ensuring a good single-cell suspension is crucial

Zorn/Bascet produces count tables and you can thus use many single-cell
tools to process the data, such as to remove doublets. Here we give an
example using **`scDblFinder`**.

## Removing doublets using scDblFinder

[scDblFinder](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html)
lives in Bioconductor and operates on a `SingleCellExperiment` object,
so the workflow is to wrap the count matrix from your Seurat object, run
the caller, then copy the per-cell scores back into the Seurat object’s
metadata.

``` r

library(SingleCellExperiment)
library(scDblFinder)

## Wrap the raw counts as a SingleCellExperiment (rows = features, cols = cells)
sce <- SingleCellExperiment::SingleCellExperiment(
  list(counts = GetAssayData(adata,layer = "counts"))
)

## Call doublets
sce <- scDblFinder::scDblFinder(sce)

## Copy the per-cell score and class back into the Seurat object
adata$scDblFinder.score <- sce$scDblFinder.score
adata$scDblFinder.class <- sce$scDblFinder.class
```

You can now overlay the score and class on the UMAP to see whether
doublets form a coherent population (often they sit between two species
clusters, or they all form a “blob” in the middle of the UMAP):

``` r

FeaturePlot(adata, features = "scDblFinder.score")
DimPlot(adata, group.by = "scDblFinder.class")
```

To remove doublets, just subset the Seurat object:

``` r

adata <- adata[, adata$scDblFinder.class == "singlet"]
```
