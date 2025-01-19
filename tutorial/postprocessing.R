
library(zorn)
library(Signac)
library(Seurat)
library(ggplot2)

### many convenience functions will be added...

################################################################################
################## Postprocessing with Signac ##################################
################################################################################


### Tell where all preprocessed data is stored, and further data will be stored
bascetRoot = "testdata"  


### Read count matrix
cnt <- ReadBascetCountMatrix(file.path(bascetRoot,"kmer.1.counts.hdf5"))  #### TODO: should read and concatenate multiple matrices; different name?

colnames(cnt) <- paste0("BASCET_",colnames(cnt)) ### compatibilíty with fragment names --- this line will likely not be needed in the future


adata <- CreateSeuratObject(
  counts = CreateKmerAssay(cnt),
  assay = "peaks"
)

DefaultAssay(adata) <- "peaks"
adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff = 'q0')
adata <- RunSVD(adata)
DepthCor(adata)
adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 2:30, reduction.name = "kmers_umap")  ## even 2nd factor is correlated

DimPlot(object = adata, label = TRUE, reduction = "kmers_umap") + NoLegend()
DimPlot(object = adata, label = TRUE, group.by = "species", reduction = "kmers_umap") + NoLegend()

FeaturePlot(adata, "nCount_peaks", reduction = "kmers_umap")




################################################################################
################## Kraken-based analysis #######################################  
################################################################################


### Read the count matrix
mat <- ReadBascetKrakenMatrix("/husky/henriksson/atrandi/wgs_miseq2/kraken_count.hdf5")


### Compress the representation to avoid trouble with some tools (will hide this code in the future)
use_row <- rowSums(mat)>0
compressed_mat <- mat[use_row, ]
rownames(compressed_mat) <- paste0("taxid-", which(use_row))


### Add count matrix to Seurat/Signac object
taxid_ob <- CreateAssayObject(compressed_mat)
adata[["kraken"]] <- taxid_ob


## Add consensus taxonomy to metadata. Will add a convenience function
kraken_taxid <- KrakenFindMaxTaxid(mat)
rownames(kraken_taxid) <- kraken_taxid$cell_id
kraken_taxid <- kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]
adata@meta.data <- cbind(adata@meta.data,kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")])



## Perform dimensional reduction using kraken count data
DefaultAssay(adata) <- "kraken"
adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff = 'q0')
adata <- RunSVD(adata)
DepthCor(adata)
adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:30, reduction.name = "kraken_umap") #1:30 or 2:30; need more testing

DimPlot(object = adata, label = TRUE, group.by = "species", reduction = "kraken_umap") + NoLegend()


## Plot species frequencies
KrakenSpeciesDistribution(adata)



################################################################################ 
################## Knee-plot per species ####################################### --- note that this is only for cells we picked!!
################################################################################ 

## based on alignment
#DefaultAssay(adata) <- "species_cnt"
#KneeplotPerSpecies(adata)

## based on kraken
DefaultAssay(adata) <- "kraken"
KneeplotPerSpecies(adata, max_species = 10)

## TODO: special call for kraken matrix, to convert to proper name? or map taxonomy ID during compression?


################################################################################
################## Barnyard plot ###############################################
################################################################################


DefaultAssay(adata) <- "species_cnt"
BarnyardPlotMatrix(adata)


## TODO: after detecting doublets, we should color by this in the plot. enable this as a feature, or plot any metadata on top



