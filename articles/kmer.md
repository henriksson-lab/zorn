# Informative KMER-based workflow

First set up Zorn/Bascet according to the [install
instructions](https://henriksson-lab.github.io/zorn/articles/install.md).
This tutorial assumes that you have [debarcoded the
reads](https://henriksson-lab.github.io/zorn/articles/debarcoding.md).

## When to use the informative KMER-based workflow

This workflow takes debarcoded reads from each cell and determines the
presence-absence of a set of KMERs for each cell. These KMERs can be
determined by a min-hash-like approach, or picked from some source of
prior knowledge. It can be used as a “second/third” workflow, after
first using the [KRAKEN2-based
workflow](https://henriksson-lab.github.io/zorn/articles/kraken.md) and
assigning taxonomic IDs, or the unbiased [Count sketch
workflow](https://henriksson-lab.github.io/zorn/articles/countsketch.md).
It can be unbiased like the Count sketch workflow, but retains
interpretability in that you will know for certain if a KMER was
detected for a cell or not. The disadvantage is that not all KMERs are
retained for the final analysis step, as it is not possible to retain
all source KMERs in memory.

## Preprocessing

To extract informative KMERs, you first need to obtain a list of them.
Bascet/Zorn offers one way, in which we first extract minhashes from
each cell. These are finger prints that we later can look for in other
cells.

[(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r

### Compute minhashes for each cell
BascetComputeMinhash(
  bascetRoot
)

### Gather minhashes into a single histogram
BascetMakeMinhashHistogram(
  bascetRoot
)
```

If two cells of the same species are shallowly sequenced, they may not
share minhashes, simply because the reads don’t span all of their
genomes. We can however look for KMERs that are more frequently shared
by counting minhashes across all cells, and generating a histogram.

``` r

kmerHist <- BascetReadMinhashHistogram(bascetRoot)

kmerHist$rank <- 1:nrow(kmerHist)
ggplot(kmerHist[sample(1:nrow(kmerHist), min(30000, nrow(kmerHist))),], aes(rank, cnt)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()
```

It is an open problem which KMERs to pick. KMERs of low abundance will
not help clustering much as the UMAP algorithm will not find many other
cells they are correlated with. Likewise, KMERs present in each cell are
not [informative
either](https://en.wikipedia.org/wiki/Binary_entropy_function). In
practice, we don’t see this happening; just picking the most common
KMERs seems to give useful results. Here we pick up to the top 100,000
KMERs and write them to a file so the query step can be rerun later:

``` r

### Pick KMERs
kmerHist$rand_index <- sample(1:nrow(kmerHist)) #add a tie breaker
kmerHist <- kmerHist[order(kmerHist$cnt, kmerHist$rand_index, decreasing = TRUE),]

useKMERs <- kmerHist$kmer[1:min(100000, nrow(kmerHist))]
writeLines(useKMERs, file.path(bascetRoot, "use_kmers.txt"))
```

Using the list of KMERs, you can now re-scan all genomes for their
presence.

[(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r

### Build count table by looking up selected KMERs in per-cell KMER databases
useKMERs <- readLines(file.path(bascetRoot, "use_kmers.txt"))
BascetQueryFq(
  bascetRoot,
  useKMERs=useKMERs
)
```

The output format is one binary HDF5 file for each shard, following the
[AnnData on-disk
layout](https://anndata.readthedocs.io/en/stable/fileformat-prose.html)
for the main sparse count matrix. Cells are stored as observations and
selected k-mers as variables. To load these, use the following command
that loads the files and concatenates them into a single matrix. It can
then be loaded into Seurat:

``` r

library(Seurat)
### Read count matrix
cnt <- ReadBascetCountMatrix(
  bascetRoot,
  "kmer_counts"
)
```

Most likely you will have to filter out the low abundant cells. Note
that with this workflow, the cutoff is based on the number of KMERs
rather than the number of reads.

``` r

  ## Load subset of real cells
  keep_cells <- rowSums(cnt)>20000 #this cutoff is based on the number of KMERs, not read count!
  sum(keep_cells)
  adata <- CreateSeuratObject(
    cnt[keep_cells,],
    assay = "infokmer"
  )
```

You can now proceed to generate a dimensional reduction in the form of a
UMAP. See the tutorials for [Seurat](https://satijalab.org/seurat/) and
[Signac](https://stuartlab.org/signac/) for details. We summarize the
required command below, where you can apply either an RNA-seq or
ATAC-seq dimensional reduction depending on what you think is
appropriate.

``` r

if(TRUE){
  #ATAC-seq style analysis
  library(Signac)
  adata <- RunTFIDF(adata)
  adata <- FindTopFeatures(adata, min.cutoff = 'q0')
  numdim <- min(30, nrow(adata)-1)
  adata <- RunSVD(adata, n = numdim)
  DepthCor(adata)
  adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:numdim, reduction.name = "infokmers_umap")
} else {
  #RNA-seq style analysis
  adata <- NormalizeData(adata)
  adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
  adata <- ScaleData(adata, features = rownames(adata))
  numdim <- min(20, nrow(adata)-2)
  adata <- RunPCA(adata, features = VariableFeatures(object = adata), npcs = numdim)
  adata <- RunUMAP(adata, dims = 1:numdim, reduction.name = "infokmers_umap")
}
```

You can now plot a UMAP with your cells. Note that the cells do not have
any labels. One way to obtain labels is by also running the KRAKEN-based
workflow; or you can use other annotation tools through the
mapcell-system.

``` r

DimPlot(object = adata, reduction = "infokmers_umap") +
  xlab("KMER1") +
  ylab("KMER2")
```
