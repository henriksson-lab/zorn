source("R/job_general.R")
source("R/job_local.R")
source("R/job_slurm.R")
source("R/bascet_file.R")
source("R/zorn.R")
source("R/shell.R")
source("R/zorn_aggr.R")
source("R/count_kmer.R")


################################################################################
################## Preprocessing with Bascet/Zorn ##############################
################################################################################

inst <- LocalInstance(direct = TRUE, show_script=TRUE)
bascetRoot = "/husky/henriksson/atrandi/wgs_miseq2/"
rawmeta <- DetectRawFileMeta("/husky/fromsequencer/240903_wgs_atcc2_miseq/raw")


### Debarcode the reads, then sort them.
BascetGetRawAtrandiWGS(
  bascetRoot,
  rawmeta,
  runner=inst
)

### Decide cells to include
h <- ReadHistogram(bascetRoot,"debarcoded")
PlotHistogram(h)
includeCells <- h$cellid[h$count>500]   ### for miseq #2
length(includeCells)

### Shardify i.e. divide into multiple sets of files for parallel processing
BascetShardify(
  bascetRoot,
  includeCells = includeCells,
  runner = inst
)

### Assemble ; 
if(FALSE){
  # assembly - 
  #test_skesa:
  #  rm -Rf temp; cargo +nightly run mapcell -i testdata/out_complete.0.tirp.gz -o testdata/skesa.0.zip -s _skesa --show-script-output
  #test_kmc_reads:
}


### Generate KMER databases for each cell
BascetMapCell(
  bascetRoot,
  withfunction = "_kmc_process_reads",
  inputName = "filtered",
  outputName = "kmc",
  runner=inst
)


### Sum up per-cell KMER databases. This will be used to produce a global histogram.
### A histogram like this can be used as /one/ way of picking KMERs of relevance
BascetFeaturise(
    bascetRoot, 
    includeCells = sample(BascetCellNames(bascetRoot, "filtered")$cell, 100),  ### avoid crash by running a subset
    runner=inst
)

### Look at the KMER histogram
KmcGetHistogram(file.path("/husky/henriksson/atrandi/wgs_miseq2","sumkmers.1.kmcdb"))  ######### todo sum up histograms

### Pick KMERs for use as features .................... some strains are very common. maybe pick across all freqs?
useKMERs <- KmcChooseKmerFeatures(
  file.path("/husky/henriksson/atrandi/wgs_miseq2","sumkmers.1.kmcdb"), 
  num_pick=10000, minfreq=0.02, maxfreq=0.30)
  


### Build count table by looking up selected KMERs in per-cell KMER databases
BascetQuery(
    bascetRoot, 
    useKMERs = useKMERs,
    runner=inst)


################################################################################
################## Postprocessing with Signac ##################################
################################################################################

### Read count matrix
cnt <- ReadBascetCountMatrix(file.path("/husky/henriksson/atrandi/wgs_miseq2","kmer.1.counts.hdf5"))  #### TODO: should read and concatenate multiple matrices; different name?

bascet_instance.default  #temp dir here


library(Signac)
library(Seurat)
library(ggplot2)


rownames(cnt) <- paste0(rownames(cnt),":1-1")

chrom_assay <- CreateChromatinAssay(
  counts = cnt, 
  sep = c(":", "-"),
#  min.cells = 10,
#  min.features = 200
)

adata <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)


adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff = 'q0')
adata <- RunSVD(adata)

DepthCor(adata)

adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 3:30)  ## even 2nd factor is correlated
#adata <- FindNeighbors(object = adata, reduction = 'lsi', dims = 2:30)
#adata <- FindClusters(object = adata, verbose = FALSE, algorithm = 3)
DimPlot(object = adata, label = TRUE) + NoLegend()

FeaturePlot(adata, "nCount_peaks")


#TODO: add metadata, depth!!!
#TODO QUAST information on top






################################################################################
################## Reference-based mapping to get "ground truth" ###############
################################################################################


