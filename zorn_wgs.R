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
source("R/kraken.R")


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



### Run Kraken on each cell  ---- these two commands should be merged
BascetRunKraken(
  bascetRoot, 
  useKrakenDB="/data/henlab/kraken/standard-8",
  numLocalThreads=10,
  runner=inst
)
BascetRunKrakenMakeMatrix(
  bascetRoot, 
  useKrakenDB="/data/henlab/kraken/standard-8",
  numLocalThreads=10,
  runner=inst
)

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

DefaultAssay(adata) <- "peaks"
adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff = 'q0')
adata <- RunSVD(adata)
DepthCor(adata)
adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 2:30, reduction.name = "kmers_umap")  ## even 2nd factor is correlated

DimPlot(object = adata, label = TRUE, reduction = "kmers_umap") + NoLegend()
DimPlot(object = adata, label = TRUE, group.by = "species", reduction = "kmers_umap") + NoLegend()

FeaturePlot(adata, "nCount_peaks", reduction = "kmers_umap")


#TODO: add metadata, depth!!!
#TODO QUAST information on top


################################################################################
################## Kraken-based analysis #######################################  
################################################################################



mat <- ReadBascetKrakenMatrix("/husky/henriksson/atrandi/wgs_miseq2/kraken_count.hdf5")


### Compress the representation to avoid trouble with some tools
use_row <- rowSums(mat)>0
compressed_mat <- mat[use_row, ]
rownames(compressed_mat) <- paste0("taxid-", which(use_row))


taxid_ob <- CreateAssayObject(compressed_mat)
adata[["kraken"]] <- taxid_ob


## Add max taxid to metadata
kraken_taxid <- KrakenFindMaxTaxid(mat)
rownames(kraken_taxid) <- kraken_taxid$cell_id
kraken_taxid <- kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]
adata@meta.data <- cbind(adata@meta.data,kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]) #AddMetaData behaves weirdly!!


#### How many species?
KrakenSpeciesDistribution(adata)


## Dimensional reduction using kraken
DefaultAssay(adata) <- "kraken"

adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff = 'q0')
adata <- RunSVD(adata)
DepthCor(adata)
adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:30, reduction.name = "kraken_umap")  ## depth may be less of a problem here; hard to tell. could normalize counts without issue. RNAseq analysis instead?
#DimPlot(object = adata, label = TRUE) + NoLegend()


DimPlot(object = adata, label = TRUE, group.by = "species", reduction = "kraken_umap")




################################################################################
################## Reference-based mapping to get "ground truth" ###############  TODO, wrap these
################################################################################


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


if(FALSE){
  BascetMapTransform(
    bascetRoot, 
    "filtered", 
    "for_bwa",
    out_format="fragments.gz",
    runner=inst)
  ### bascet bam2fragments 
}

################################################################################
################## Reference-based mapping to get "ground truth" ###############
################################################################################


### Get reads in fastq format  --- will not be needed later
BascetMapTransform(
  bascetRoot, 
  "filtered", 
  "asfq",
  out_format="fq.gz",   ## but we need two fq as out!! ideally at least. or if R1.fq.gz => write two of them. otherwise gather?
  runner=inst
)

### Perform alignment -- this is a wrapper for mapshard
BascetAlignToReference(
  bascetRoot,
  useReference="/husky/fromsequencer/241016_novaseq_wgs2/trimmed/ref10/all.fa",
  numLocalThreads=10
)

### Generate fragments BED file suited for quantifying reads/chromosome using Signac later -- this is a wrapper for mapshard
BascetBam2Fragments(
  bascetRoot,
  runner=inst
)

  
################################################################################
################## Get alignment stats per cell ################################   
################################################################################


## Add the counts to our Seurat object
adata[["chrom_cnt"]] <- FragmentCountsPerChromAssay(bascetRoot)
DefaultAssay(adata) <- "chrom_cnt"

#Figure out which chromosome has most reads in which cell
cnt <- adata@assays$chrom_cnt$counts
adata$dominant_chr <- rownames(cnt)[apply(cnt, 2, which.max)]

DimPlot(adata, group.by = "dominant_chr")


# #NC_017316.1   Enterococcus faecalis OG1RF, complete sequence
# #CP086328.1    Bacillus pacificus strain anQ-h4 chromosome, complete genome



################################################################################ 
################## Which sequences belong to the same strain? ##################
################################################################################ 

map_strain2gram <- read.csv("/husky/fromsequencer/240701_wgs_atcc1/straintype.csv")

library(stringr)

get_map_seq2strain <- function(){
  map_seq2strain <- NULL
  refdir <- "/husky/fromsequencer/240809_novaseq_wgs1/trimmed/ref10/separate"
  for(f in list.files(refdir, pattern = "*.fasta")){
    print(f)
    thel <- readLines(file.path(refdir,f))
    thel <- thel[str_starts(thel,">")]
    onedf <- data.frame(line=thel)
    onedf$strain <- f
    map_seq2strain <- rbind(map_seq2strain, onedf)
  }
  map_seq2strain$id <- str_split_fixed(str_sub(map_seq2strain$line,2)," ",2)[,1]
  map_seq2strain$name <- str_split_fixed(str_sub(map_seq2strain$line,2)," ",2)[,2]
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain,".fasta")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," chr1")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," chr2")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," pRSPH01")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," pMP1")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," pCP1")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," NC_003909.8")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," NC_005707.1")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," NC_005707.1")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," plasmid")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," CP086328.1")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," CP086329.1")

  map_seq2strain <- merge(map_seq2strain,map_strain2gram)
  map_seq2strain <- map_seq2strain[,colnames(map_seq2strain)!="line"]
  
  idx <- read.table(pipe("samtools idxstats /husky/fromsequencer/241206_novaseq_wgs3/trimmed/sorted.bam"),sep="\t")
  colnames(idx) <- c("id","len","mapped","unmapped")
  map_seq2strain <- merge(idx[,c("id","len")], map_seq2strain)  #note, only using id and len
  
  map_seq2strain
}




map_seq2strain <- get_map_seq2strain()[,c("id","len","strain")]
strain_genomesize <- sqldf::sqldf("select sum(len) as len, strain from map_seq2strain group by strain")


########## Produce a count matrix on strain level
DefaultAssay(adata) <- "chrom_cnt"
adata[["species_cnt"]] <- ChromToSpeciesCount(adata, map_seq2strain)

#Figure out which species has most reads in which cell
cnt <- adata@assays$species_cnt$counts
adata$dominant_species <- rownames(cnt)[apply(cnt, 2, which.max)]

DimPlot(adata, group.by = "species_aln")





################################################################################ 
################## Knee-plot per species ####################################### --- note that this is only for cells we picked!!
################################################################################ 


DefaultAssay(adata) <- "species_cnt"
KneeplotPerSpecies(adata)

DefaultAssay(adata) <- "kraken"
KneeplotPerSpecies(adata, max_species = 10)

## TODO: special call for kraken matrix, to convert to proper name?  ; map taxid during compression!





################################################################################
################## Barnyard plot ###############################################
################################################################################



DefaultAssay(adata) <- "species_cnt"
BarnyardPlotMatrix(adata)


## TODO: after detecting doublets, we should color by this in the plot   ; color by any metadata!!





################################################################################
################## Correlation of genomes ######################################
################################################################################

DefaultAssay(adata) <- "species_cnt"


SpeciesCorrMatrix(adata)


# other code in /home/mahogny/jupyter/scWGS







