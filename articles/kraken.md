# KRAKEN2

First set up Zorn/Bascet according to the [install
instructions](https://henriksson-lab.github.io/zorn/articles/install.md).
This tutorial assumes that you have [debarcoded the
reads](https://henriksson-lab.github.io/zorn/articles/debarcoding.md).

## When to use the KRAKEN2 workflow

This workflow takes debarcoded reads and assigns taxonomic IDs to each
of them. This results in a count matrix that can be used for dimensional
reduction and clustering. The advantage of this workflow is that (1) it
lets you use prior information to group cells, (2) it lets you assign a
taxonomic ID (species) to each cell using majority voting. The workflow
is fast and a good “first QC” for any metagenomic dataset

## Classifying reads with KRAKEN2

To run KRAKEN2, you need a database. You can get them here:
<https://benlangmead.github.io/aws-indexes/k2>

If you want a small one then consider standard-8. Unzip it in a
directory. You can then run KRAKEN2 like this:

[(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r

### Run Kraken on each cell. Produce a count matrix of taxonomy features
BascetRunKraken(
  bascetRoot,
  useKrakenDB="/your_disk/kraken/standard-8"
)
```

Internally, there are two steps here. First KRAKEN2 taxonomically
classifies each read. In the second step, we count the taxonomic reads
for each cell. This ends up being a rather small count matrix that you
can process using Seurat.

## Postprocessing with Signac/Seurat

The output format is one binary HDF5 file for each shard, following the
[AnnData on-disk
layout](https://anndata.readthedocs.io/en/stable/fileformat-prose.html)
for the main sparse count matrix. Cells are stored as observations,
taxonomy IDs as variables, and per-cell unclassified read counts are
available in `mat@obs$unclassified_reads`. To load these, use the
following command that loads the files and concatenates them into a
single matrix. It can then be loaded into Seurat:

``` r

mat <- ReadBascetCountMatrix(bascetRoot,"kraken_mat", verbose=FALSE)

adata <- CreateSeuratObject(mat, project = "proj", min.cells = 0, min.features = 0)
```

You then most likely want to know which cell is which species. Zorn can
annotate species by simply picking the dominant species, phylum and
group for each cell:

``` r

## Look up taxonomy consensus data given taxonomyID counts for each cell
kraken_taxid <- KrakenFindConsensusTaxonomy(mat)

## Add KRAKEN consensus taxonomy to metadata
adata <- BascetAddMetaData(adata, kraken_taxid)
```

You will need to filter low-abundance cells. To do this, first
investigate a kneeplot of each species. As different species are
differently hard to lyse, you will likely have phylum-specific kneeplot
patterns.

``` r

showNumSpecies <- 10
KrakenKneePlot(adata, groupby = "phylum", showNumSpecies=showNumSpecies)
```

You can then proceed to filter out cells with few reads. Note that a
general cutoff is nonideal as different phyla have different amount of
DNA (due to different genome sizes & lysis biases); the best way to
handle this is an open research question (e.g. you could also consider
different cutoffs for different species). This is however the most basic
filtering you can perform:

``` r

keep_cells <- adata$nCount_RNA > 10000 #10k reads
sum(keep_cells)
adata <- adata[,keep_cells]
```

You can now proceed to generate a dimensional reduction in the form of a
UMAP. See the tutorials for [Seurat](https://satijalab.org/seurat/) and
[Signac](https://stuartlab.org/signac/) for details. We summarize the
required command below, where you can apply either an RNA-seq or
ATAC-seq dimensional reduction depending on what you think is
appropriate.

ATAC-seq style analysis

``` r

#ATAC-seq style analysis
library(Signac)
adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff = 'q0')
adata <- RunSVD(adata)
DepthCor(adata)
adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:30, reduction.name = "kraken_umap")  
```

RNA-seq style analysis

``` r
#RNA-seq style analysis
adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
adata <- ScaleData(adata, features = rownames(adata))
adata <- RunPCA(adata, features = VariableFeatures(object = adata))
adata <- RunUMAP(adata, dims = 1:20, reduction.name = "kraken_umap")
}
```

You can now plot a UMAP with your cells! Here colored by genus, but you
can also try other levels of annotation:

``` r

DimPlot(object = adata, label = TRUE, group.by = "genus", reduction = "kraken_umap") + 
  NoLegend() + 
  xlab("KRAKEN1") + 
  ylab("KRAKEN2")
```
