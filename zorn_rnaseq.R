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
if(TRUE){
  bascetRoot = "/husky/henriksson/atrandi/rnaseq3/1/"
  rawmeta <- DetectRawFileMeta("/husky/fromsequencer/250108_joram_rnaseq3/raw/miseq_demul/1/")
} else {
  bascetRoot = "/husky/henriksson/atrandi/rnaseq3/2/"
  rawmeta <- DetectRawFileMeta("/husky/fromsequencer/250108_joram_rnaseq3/raw/miseq_demul/2/")
}



bascet_instance.default  #temp dir is here


### Debarcode the reads, then sort them.
BascetGetRaw(
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
  #includeCells = includeCells, ## include all cells! filter later   TODO it is really expensive to search for individual cells. need other mode for merging
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
  useReference="/husky/fromsequencer/241210_joram_rnaseq/ref/all.fa",
  numLocalThreads=10
)

#"/husky/fromsequencer/241210_joram_rnaseq/ref/all.fa"

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
gff <- gff[gff$type != "exon",]
gff <- gff[gff$type != "CDS",]
gff <- gff[gff$type != "region",]


### all ID??
grange_gene <- GenomicRanges::makeGRangesFromDataFrame(gff)
grange_gene$type <- gff$type
grange_gene$gene_biotype <- gff$gene_biotype
grange_gene$Name <- gff$Name
grange_gene$Name[is.na(grange_gene$Name)] <- gff$ID[is.na(grange_gene$Name)]



#bascetRoot
#"/husky/henriksson/atrandi/rnaseq3/1/fragments.1.tsv.gz"

####### Perform counting
adata_f <- FragmentsToSignac(file.path(bascetRoot,"fragments.1.tsv.gz"))  ## TODO, support multiple fragment files
#adata[["RNA"]] <- CountGrangeFeatures(adata, grange_gene) ## this fails for no good reason; replacement has 8799 rows, data has 2
adata <- CreateSeuratObject(CountGrangeFeatures(adata_f, grange_gene))





####### 
####### Divide libraries
####### 

adata$well_r1 <- stringr::str_split_i(colnames(adata),"_",2)

AtoH <- c("A","B","C","D","E","F","G","H")
libname <- c(rep("lib1",8),rep("lib2",8),rep("lib3",8))
names(libname) <- c(paste0(AtoH,1), paste0(AtoH,2), paste0(AtoH,3))
adata$libname <- unname(libname[adata$well_r1])


####### 
####### Try to see which SPCs have cells
####### 

if(FALSE){
  #https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html#detecting-empty-droplets
  BiocManager::install("DropletUtils")
}

#my.counts should be sparse matrix before filtering

#colSums(adata@assays$RNA@counts)

e.out <- DropletUtils::emptyDrops(adata@assays$RNA@counts, lower = 50)  #how to think about this param?
adata$is_cell <- e.out$FDR <= 0.01
adata$is_cell_logprob <- e.out$LogProb
adata$is_cell_cnt <- e.out$Total
ggplot(adata@meta.data, aes(is_cell_cnt, -is_cell_logprob, color=is_cell)) + geom_point()



####### 
####### Perform knee plot
####### 

### All cells together
df <- data.frame(cnt=sort(adata$nCount_RNA, decreasing = TRUE))
df$rank <- 1:nrow(df)
ggplot(df,aes(rank, cnt)) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10()

### for each library
all_df <- NULL
for(curlib in unique(adata$libname)){
  df <- data.frame(cnt=sort(adata$nCount_RNA[adata$libname==curlib], decreasing = TRUE), lib=curlib)
  df$rank <- 1:nrow(df)
  all_df <- rbind(all_df, df)
}
ggplot(all_df,aes(rank, cnt, color=lib)) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10()



####### 
####### Perform RNA-seq type analysis; see https://satijalab.org/seurat/articles/adata3k_tutorial.html
####### 

DefaultAssay(adata) <- "RNA"

#hist(log10(adata$nCount_RNA))
adata <- adata[,adata$nCount_RNA>10] #meaningless below for certain!!

adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(adata), 10)

#[1] "npr"              "glmM"             "rpoE"             "tsf"              "rplB"             "hupB"            
#[7] "bamA"             "lpxC"             "ssrA"             "rna-BZ22-RS22510"

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
adata <- RunUMAP(adata, dims = 1:10)

adata <- FindClusters(adata, resolution = 0.2)

DimPlot(adata, reduction = "umap")

FeaturePlot(adata, features=top10)





# library(dplyr)
# pbmc.markers <- FindAllMarkers(adata, only.pos = FALSE)
# 
# df <- as.data.frame(pbmc.markers %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 0.2) %>%
#   dplyr::filter(pct.1 > 0.01))
# df$Name <- df$gene
# df <- merge(df, gff[,c("Name","type")])
# df <- df[order(df$p_val),]
# df
# 
# 
# pbmc.markers <- FindMarkers(adata, only.pos = FALSE, ident.1=0, ident.2=1, min.pct = 0.02)
# pbmc.markers
# 
# 
# 
# print(n=100,pbmc.markers %>%
#         group_by(cluster) %>%
#         dplyr::filter(pct.1 > 0.1))


FeaturePlot(adata, features=c("ssrA"))
FeaturePlot(adata, features=c("yqhD"))
FeaturePlot(adata, features=c("nlpD"))
FeaturePlot(adata, features=c("ompA"))
FeaturePlot(adata, features=c("phoH"))


FeaturePlot(adata, features=c("katG"))
FeaturePlot(adata, features=c("ssrA")) #drives a large part of the clustering. transfer-messenger RNA"
FeaturePlot(adata, features=c("BZ22-RS01385")) #another major driver. rRNA

FeaturePlot(adata, features=c("BZ22-RS04865"))



####### 
####### QC for each library (1-3)
####### 

DimPlot(adata, reduction = "umap", group.by = "libname")  #lib2 enriched in ssrA cluster

DimPlot(adata, reduction = "umap", group.by = "well_r1")#, label = TRUE)



####### 
####### For each library, count features
####### 

#curlib <- "lib1"
alldf <- NULL
for(curlib in paste0("lib",1:3)){
  df <- as.data.frame(merge(
    data.frame(
      Name=rownames(adata),
      cnt=rowSums(adata[,adata$libname==curlib]@assays$RNA@counts),
      lib=curlib
    ),
    as.data.frame(grange_gene)[,c("Name","gene_biotype")]
  ))
  #as.data.frame(df)
  alldf <- rbind(
    alldf,
    sqldf::sqldf("select gene_biotype, sum(cnt) as cnt, lib from df group by gene_biotype")
  )
  #protein_coding is the relevant one  
}
sum_cnt_type <- reshape2::acast(data = alldf, gene_biotype~lib, value.var = "cnt")

sum_cnt_type_norm <- sum_cnt_type
for(i in 1:3){
  sum_cnt_type_norm[,i] <- round(digits = 3,sum_cnt_type_norm[,i]/sum(sum_cnt_type_norm[,i]))
}

sum_cnt_type
sum_cnt_type_norm




####### 
####### Pileups
####### 

adata_f <- CreateSeuratObject(adata_f)
adata_f

adata_f$nCount_RNA <- adata$nCount_RNA
adata_f$libname <- adata$libname

CoveragePlot(
#  assay.scale = "separate",
  object = adata_f,
  region = "NZ_CP009792.1-60000-90000", 
  group.by = "libname"
)


