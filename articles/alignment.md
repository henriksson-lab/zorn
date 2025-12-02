# Alignment-based workflow

First set up your Zorn/Bascet workdirectory as before. If you wish to
run these steps on a SLURM cluster, see separate vignette and adapt
accordingly.

``` r
library(Zorn)
bascet_runner.default <- LocalRunner(direct = TRUE, showScript=TRUE)
bascetRoot <- "/home/yours/an_empty_workdirectory"
```

## Alignment

We will use BWA for alignment, which wants the data in FASTQ format.
Here we assume that you use a single-cell WGS chemistry that produces
paired reads, and direct Zorn to produce R1/R2 FASTQ files (only R1
specified). We name the output Bascets “asfq”, taking as input the
typical sharded reads:

[(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r
### Get reads in fastq format
BascetMapTransform(
  bascetRoot,
  inputName="filtered",   #default; can omit
  outputName="asfq",  ###name parameter
  out_format="R1.fq.gz"
)
```

Before alignment, you need to prepare your reference by indexing it. You
can do this using “bwa index” in the terminal; but you can also do it
via Zorn. This assumes that you have a FASTA-file named all.fa (can also
be gzipped):

[(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r
#Build a reference
BascetIndexGenomeBWA(
  bascetRoot,
  "/your_reference/all.fa"
)
```

You can now proceed with alignment using BWA. By default, this command
outputs both unsorted and sorted BAM-files.

- The sorted BAM-files are required for counting of reads in genomic
  regions, i.e., WGS species counting and RNA-seq
- The unsorted BAM-files are required for filtering of host reads

[(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r
### Perform alignment
BascetAlignToReference(
  bascetRoot,
  useReference="/your_reference/all.fa",
  numLocalThreads=10
)
```

TODO: Filtering

## Host filtering workflow

``` r
BascetFilterAlignment(
  bascetRoot,
  outputName="filter_paired",
  keep_mapped=TRUE
)
```

## Option \#1: Counting with Signac

If you want to analyze where reads mapped, either for RNA-seq counting
for genes, or more generally, what the count is within a genomic region,
you can use Signac for great flexibility in this task. Signac is
designed to analyze single-cell ATAC-seq data, where the reads are
stored in a specialized BED-file. You can use the command below to
generate this file.

Note one caveat: depending on chemistry, this workflow may not handle
UMIs, and just report raw read counts. However, for quick and dirty
operations this is typically enough.

[(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r
### Generate fragments BED file suited for quantifying
### reads/chromosome using Signac later 
BascetBam2Fragments(
  bascetRoot
)
```

## Option \#2: Counting chromosomes (species)

A simplified scenario is when you just want to get the number of reads
per chromosome. These counts can later be summarized on the levels of
genomes (summing up across chromosomes). We use this to for example
compare the alignment to an expected mock community, and thus just want
to see how species mix together. The following command generates counts
files:

[(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r
### Count reads per chromosome
BascetCountChrom(
  bascetRoot
)
```

The output format is one binary HDF5 file for each shard, roughly in the
[Anndata file
format](https://anndata.readthedocs.io/en/stable/fileformat-prose.html).
To load these, use the following command that loads the files and
concatenates them into a single matrix. It can then be loaded into
Seurat:

``` r
library(Seurat)

cnt <- ReadBascetCountMatrix(bascetRoot,"chromcount", verbose=FALSE)

adata <- CreateSeuratObject(
  counts = CreateAssayObject(t(cnt$X)), ## according to anndata standard, there should be a transpose here
  assay = "chrom_cnt"
) 
```

Most likely, some cells will have to be discarded. Initially we
recommend that you keep the cutoff low until you are better informed.
You can then apply the cutoff as follows:

``` r
keep_cells <- adata$nCount_chrom_cnt > 10000 #10k reads
sum(keep_cells)                              #See how many cells pass this cutoff
adata <- adata[,keep_cells]                  #Reduce to sensible number
```

If you have a list of which chromosomes belong together, then you can
furthermore sum them up to get species-level counts. mapSeq2strain
should contain two columns: (todo) (TODO filter after instead)

``` r
adata[["species_cnt"]] <- ChromToSpeciesCount(adata, mapSeq2strain)  

#Optional: Figure out which species has most reads in which cell
cnt <- adata@assays$species_cnt$counts
adata$species_aln <- rownames(cnt)[apply(cnt, 2, which.max)]
```

Knowing the species, you can generate a per-species kneeplot as follows:

``` r
DefaultAssay(adata) <- "species_cnt"
KneeplotPerSpecies(adata)
```

You can now proceed to generate a dimensional reduction in the form of a
UMAP. See the tutorials for [Seurat](https://satijalab.org/seurat/) and
[Signac](https://stuartlab.org/signac/) for details. We summarize the
required command below, where you can apply either an RNA-seq or
ATAC-seq dimensional reduction depending on what you think is
appropriate.

``` r
if(FALSE){
  #ATAC-seq style analysis
  library(Signac)
  adata <- RunTFIDF(adata)
  adata <- FindTopFeatures(adata, min.cutoff = 'q0')
  adata <- RunSVD(adata)
  DepthCor(adata)
  adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:(nrow(adata)-1), reduction.name = "chrom_umap")  
} else {
  #RNA-seq style analysis
  adata <- NormalizeData(adata)
  adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = nrow(adata))
  adata <- ScaleData(adata, features = rownames(adata))
  adata <- RunPCA(adata, features = VariableFeatures(object = adata))
  adata <- RunUMAP(adata, dims = 1:(nrow(adata)-1), reduction.name = "chrom_umap")
}
```

Finally, you can plot your UMAP. If you also extracted the dominant
species column in your matrix, you can visualize it on top!

``` r
DimPlot(object = adata, label = TRUE, group.by = "species_aln") + 
  xlab("BWA1") + 
  ylab("BWA2")
```

## Option \#3: RNA-seq analysis (feature counting)

For protocols such as RNA-seq, the aim is to get the number of reads per
gene, per cell. After alignment and sorting of reads based on genome
coordinate (default setting), you can overlay the reads with gene
coordinates from a GTF or GFF file (both are supported). The following
command performs this operation:

``` r
BascetCountFeature(
  bascetRoot,
  gffFile = "/my/genome.gtf.gz",
  useFeature = "gene",     #Default
  attrGeneId = "gene_id",  #Default
  attrGeneName = "name"    #Default
)
```

The ID and name of genes are extracted from the GFF attributes column.
Not all files have both name and ID; the ID will be used as name if
nothing specified. By default, the rows having “gene” as the feature
description will

After the features (genes) have been counted, you can load them into
Seurat:

``` r
library(Seurat)

cnt <- ReadBascetCountMatrix(bascetRoot,"featurecount", verbose=FALSE)

adata <- CreateSeuratObject(
  counts = CreateAssayObject(t(cnt$X)) ## according to anndata standard, there should be a transpose here
) 
```

From here, you can follow an [RNA-seq tutorial for
Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial). But for
the impatient, this is some minimal example code to get your first UMAP:

``` r
adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata, selection.method = "vst")
adata <- ScaleData(adata)
adata <- RunPCA(adata, features = VariableFeatures(object = adata))

adata <- FindNeighbors(adata, dims = 1:20)
adata <- FindClusters(adata, resolution = 0.5)

adata <- RunUMAP(adata, dims=1:10)
```
