library(zorn)

source("R/job_general.R")
source("R/job_local.R")
source("R/job_slurm.R")
source("R/bascet_file.R")
source("R/zorn.R")
source("R/shell.R")
source("R/zorn_aggr.R")
source("R/count_kmer.R")
source("R/refgenome.R")


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



### Generate min-hash from each per-cell KMER database
BascetMapCell(
  bascetRoot,
  withfunction = "_minhash",
  inputName = "kmc",
  outputName = "minhash",
  runner=inst
)





if(TRUE){
  ############################ Could be good to use a #read cutoff; this avoids pulling in noisy kmers. how about weight 1/#read count?
  all_kmer <- AggregateMinhashes(bascetRoot) 
  
  PlotMinhashDistribution <- function(all_kmer){
    ggplot(all_kmer, aes(index, freq)) + 
      geom_line() + 
      scale_x_log10() + 
      scale_y_log10() + theme_bw()
  } 
  PlotMinhashDistribution(all_kmer)

  
  useKMERs <- all_kmer$kmer[all_kmer$freq>3]
  length(useKMERs)
  
}




if(FALSE){
  
  
  ### Sum up per-cell KMER databases. This will be used to produce a global histogram. ---------------- this function can be taken out
  ### A histogram like this can be used as /one/ way of picking KMERs of relevance
  BascetFeaturise(
    bascetRoot, 
    includeCells = sample(BascetCellNames(bascetRoot, "filtered")$cell, 100),  ### avoid crash by running a subset
    runner=inst
  )
  
  ### Look at the KMER histogram
  KmcGetHistogram(file.path(bascetRoot,"sumkmers.1.kmcdb"))  ######### todo sum up histograms
  
  ### Pick KMERs for use as features .................... some strains are very common. maybe pick across all freqs?
  useKMERs <- KmcChooseKmerFeatures(
    file.path(bascetRoot,"sumkmers.1.kmcdb"), 
    num_pick=10000, minfreq=0.02, maxfreq=0.30)  
}

  


### Build count table by looking up selected KMERs in per-cell KMER databases
BascetQuery(
    bascetRoot, 
    useKMERs = useKMERs,
    runner=inst)




################################################################################
################## Postprocessing with Signac ##################################
################################################################################

library(Signac)
library(Seurat)
library(ggplot2)


### Read count matrix
cnt <- ReadBascetCountMatrix(file.path(bascetRoot,"kmer.1.counts.hdf5"))  #### TODO: should read and concatenate multiple matrices; different name?

#bascet_instance.default  #temp dir here

colnames(cnt) <- paste0("BASCET_",colnames(cnt)) ### compatibilíty with fragments


adata <- CreateSeuratObject(
  counts = CreateKmerAssay(cnt),
  assay = "peaks"
)


adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff = 'q0')
adata <- RunSVD(adata)

DepthCor(adata)

adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 2:30)  ## even 2nd factor is correlated
DimPlot(object = adata, label = TRUE) + NoLegend()

FeaturePlot(adata, "nCount_peaks")


#TODO: add metadata, depth!!!
#TODO QUAST information on top






################################################################################
################## Reference-based mapping to get "ground truth" ###############
################################################################################


### Get reads in fastq format  --- will not be needed later
BascetMapTransform(
  bascetRoot, 
  "filtered", 
  "asfq",
  out_format="fq.gz",   ## but we need two fq as out!! ideally at least. or if R1.fq.gz => write two of them. otherwise gather?
  runner=inst)


### Perform alignment via mapshard
### not yet implemented. should be done in bascet, not zorn, to support piping
if(FALSE){
  BascetMapShard(
    bascetRoot,
    withfunction = "_bwa",
    inputName = "for_bwa",
    outputName = "aligned",
    runner=inst
  )
}

### Perform alignment -- this is a wrapper for mapshard
BascetAlignToReference(
  bascetRoot,
  useReference="/husky/fromsequencer/241016_novaseq_wgs2/trimmed/ref10/all.fa",
  numLocalThreads=10
)

### Generate BED file suited for quantifying features using Signac later -- this is a wrapper for mapshard
BascetBam2Fragments(
  bascetRoot,
  runner=inst
)


if(FALSE){
  BascetMapTransform(
    bascetRoot, 
    "filtered", 
    "for_bwa",
    out_format="fragments.gz",
    runner=inst)
  ### bascet bam2fragments 
}
  
  
### Now possible to get read count, per cell, for each chromosome!!
  
  

################################################################################
################## Get alignment stats per cell ################################   
################################################################################


## Add the counts to our Seurat object
adata[["chrom_cnt"]] <- FragmentCountsPerChromAssay(bascetRoot)

DefaultAssay(adata) <- "chrom_cnt"

cnt <- adata@assays$chrom_cnt$counts
adata$dominant_chr <- rownames(cnt)[apply(cnt, 2, which.max)]


DimPlot(adata, group.by = "dominant_chr")


# #NC_017316.1   Enterococcus faecalis OG1RF, complete sequence
# #CP086328.1    Bacillus pacificus strain anQ-h4 chromosome, complete genome



