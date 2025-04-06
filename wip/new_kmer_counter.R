
bascet_runner <- LocalRunner(direct = TRUE, show_script=TRUE)
bascetRoot = "/husky/henriksson/atrandi/v2_wgs_miseq/"
#bascetRoot = "/husky/henriksson/atrandi/wgs_novaseq3/"
#rawmeta <- DetectRawFileMeta("/husky/fromsequencer/240903_wgs_atcc2_miseq/raw")


### Compute minhashes for each cell
BascetComputeMinhash(
  bascetRoot,
  runner=bascet_runner
)

### Gather minhashes into a single histogram
BascetMakeMinhashHistogram(
  bascetRoot,
  runner=bascet_runner
)

### Pick KMERs
kmer_hist <- BascetReadMinhashHistogram(bascetRoot)
useKMERs <- kmer_hist$kmer[kmer_hist$cnt>5]

kmer_hist$rank <- 1:nrow(kmer_hist)
ggplot(kmer_hist[kmer_hist$cnt>2,], aes(rank, cnt)) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10()


### Build count table by looking up selected KMERs in per-cell KMER databases
BascetQueryFq(
  bascetRoot,
  useKMERs=useKMERs,
  runner=bascet_runner
)

################################################################################
################## Postprocessing with Signac ##################################
################################################################################

library(Signac)
library(Seurat)
library(ggplot2)


### Read count matrix
cnt <- ReadBascetCountMatrix(file.path(bascetRoot,"kmer_counts.h5"))  #### TODO: should read and concatenate multiple matrices; different name?
#cnt <- ReadBascetCountMatrix(file.path(bascetRoot,"kmer_5.h5"))  #### TODO: should read and concatenate multiple matrices; different name?


adata <- CreateSeuratObject(
  counts = CreateAssayObject(cnt),
  assay = "peaks"
)

DefaultAssay(adata) <- "peaks"
adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff = 'q0')
adata <- RunSVD(adata)
DepthCor(adata)
adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 2:30, reduction.name = "kmers_umap")  ## miseq, 1st factor is ok!ish

FeaturePlot(adata, "nCount_peaks", reduction = "kmers_umap")




