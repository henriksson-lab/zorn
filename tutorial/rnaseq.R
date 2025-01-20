library(Zorn)

if(FALSE){
  source("R/job_general.R")
  source("R/job_local.R")
  source("R/job_slurm.R")
  source("R/bascet_file.R")
  source("R/zorn.R")
  source("R/shell.R") 
  source("R/zorn_aggr.R")
  source("R/count_kmer.R")
  source("R/refgenome.R")
}


################################################################################
################## Preprocessing with Bascet/Zorn ##############################
################################################################################

inst <- LocalInstance(direct = TRUE, show_script=TRUE)
bascetRoot = "/husky/henriksson/atrandi/rnaseq3/1/"
rawmeta <- DetectRawFileMeta("/husky/fromsequencer/250108_joram_rnaseq3/raw/miseq_demul/1/")

bascet_instance.default  #temp dir is here


### Debarcode the reads, then sort them.
BascetGetRawAtrandiWGS(
  bascetRoot,
  rawmeta,
  runner=inst
)

### Decide cells to include
h <- ReadHistogram(bascetRoot,"debarcoded")
PlotHistogram(h)
includeCells <- h$cellid[h$count>100]   ### for rnaseq miseq #3.1
length(includeCells)

### Shardify i.e. divide into multiple sets of files for parallel processing
BascetShardify(
  bascetRoot,
  includeCells = includeCells,
  runner = inst
)



### Note: shardification can partially be skipped for RNAseq if you know what you are doing.
### but this is poorly tested, so doing it to be sure!






################################################################################
################## Reference-based mapping to get "ground truth" ###############
################################################################################


### Get reads in fastq format for BWA  --- will not be needed later as the conversion will be made on the fly
BascetMapTransform(
  bascetRoot, 
  "filtered", 
  "asfq",
  out_format="fq.gz",   ## TODO we need two fq as out!! ideally at least. or if R1.fq.gz => write two of them. otherwise gather?
  runner=inst
)

### Perform alignment -- internally wraps mapshard
BascetAlignToReference(
  bascetRoot,
  useReference="/husky/fromsequencer/241016_novaseq_wgs2/trimmed/ref10/all.fa",
  numLocalThreads=10
)

### Generate fragments BED file suited for quantifying reads/chromosome using Signac later
BascetBam2Fragments(
  bascetRoot,
  runner=inst
)


################################################################################
################ Use Signac/Seurat as a feature counter ########################
################################################################################



library(Signac)
library(Seurat)
library(ggplot2)


####### Find coordinates of each gene for our organism of choice
gff <- rtracklayer::readGFF("/husky/fromsequencer/241210_joram_rnaseq/ref/all.gff3")
gff_gene <- gff[gff$type=="gene",]
gff_gene <- unique(gff_gene[,c("seqid","start","end","Name","gene_biotype")])
grange_gene <- GenomicRanges::makeGRangesFromDataFrame(gff_gene)
grange_gene$Name <- gff_gene$Name
grange_gene$gene_biotype <- str_replace_all(gff_gene$gene_biotype,"_","-")


####### Perform counting
adata <- LoadUncountedFragments("/husky/fromsequencer/241210_joram_rnaseq/trimmed/atac_fragments.tsv.gz")  ## TODO, support multiple fragment files
adata[["RNA"]] <- CountGrangeFeatures(grange_gene)


####### 
####### Perform RNA-seq type analysis; see https://satijalab.org/seurat/articles/adata3k_tutorial.html
####### 

DefaultAssay(adata) <- "RNA"


adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(adata), 10)

# plot variable features with and without labels
if(FALSE){
  plot1 <- VariableFeaturePlot(adata)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
}


all.genes <- rownames(adata)
adata <- ScaleData(adata, features = all.genes)
adata <- RunPCA(adata, features = VariableFeatures(object = adata))
adata <- FindNeighbors(adata, dims = 1:10)
adata <- FindClusters(adata, resolution = 0.5)
adata <- RunUMAP(adata, dims = 1:10)

DimPlot(adata, reduction = "umap")





#TODO add rRNA info!


