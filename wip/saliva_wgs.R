
################################################################################
################## Preprocessing with Bascet/Zorn ##############################
################################################################################

#bascet_inst <- LocalInstance(direct = TRUE, show_script=TRUE)
bascet_runner.default <- LocalRunner(direct = TRUE, show_script=FALSE)
bascet_instance.default <- BascetInstance(bin = "/home/mahogny/github/bascet/target/debug/bascet", tempdir = "/tmp/")
bascet_instance.default <- getBascetSingularityImage("/home/mahogny/github/bascet/singularity/") #BascetInstance(bin = "/home/mahogny/github/bascet/target/debug/bascet", tempdir = "/tmp/")#, direct = TRUE, show_script=FALSE)

bascetRoot <- "/husky/henriksson/atrandi/wgs_saliva1/"
rawmeta <- DetectRawFileMeta("/husky/fromsequencer/250311_failed_saliva/saliva_wgs/raw")


### Debarcode the reads, then sort them.  BUG, BGZIP IS ONLY USING ONE CPU
BascetGetRaw(
  bascetRoot,
  rawmeta
)

### Decide cells to include
h <- ReadHistogram(bascetRoot,"debarcoded")
PlotHistogram(h)
includeCells <- h$cellid[h$count>500]   ### for miseq #2
length(includeCells)

### Shardify i.e. divide into multiple sets of files for parallel processing  BUG ONLY USING ONE CPU
BascetShardify(
  bascetRoot
  #  includeCells = includeCells,  #TODO: best to NOT filter anything at this stage. use includeCells in later stage for quality genomes
  #  runner = bascet_inst
)

################################################################################
################## Preprocessing with KRAKEN ###################################
################################################################################

## By design, running KRAKEN on the full list of reads is pretty cheap. Doing so enables
## automatic empty droplet calling using DropletUtils later


### Get reads in fastq format  --- will not be needed later
BascetMapTransform(
  bascetRoot, 
  "filtered", 
  "asfq",
  out_format="fq.gz"   ## but we need two fq as out!! ideally at least. or if R1.fq.gz => write two of them. otherwise gather?
  #  runner=bascet_inst
)


### Run Kraken on each cell  ---- these two commands should be merged
BascetRunKraken(  ### BUG!!! detected 4 input files, but only processed the first
  bascetRoot, 
  useKrakenDB="/data/henlab/kraken/standard-8",
  numLocalThreads=10
  #  runner=bascet_inst
)
BascetRunKrakenMakeMatrix(
  bascetRoot, 
  useKrakenDB="/data/henlab/kraken/standard-8",
  numLocalThreads=10
  #  runner=bascet_inst
)






################################################################################
################## Kraken-based analysis #######################################  
################################################################################

#bascetRoot <- "/husky/henriksson/atrandi/wgs_novaseq3/"




### TODO -- need to assemble 4 kraken matrices

mat <- ReadBascetKrakenMatrix(file.path(bascetRoot,"kraken.1.counts.hdf5")) 
## any reason this need to be separate from the usual count table reading? best if not!


### Compress the representation to avoid trouble with some tools
compressed_mat <- SetTaxonomyNamesFeatures(mat)

taxid_ob <- CreateAssayObject(compressed_mat)
adata <- CreateSeuratObject(counts = taxid_ob, project = "proj", min.cells = 0, min.features = 0) ### do we need this? can overload also on assayobject!

#adata[["kraken"]] <- taxid_ob


## Add KRAKEN consensus taxonomy to metadata
kraken_taxid <- KrakenFindConsensusTaxonomy(mat)
rownames(kraken_taxid) <- kraken_taxid$cell_id
kraken_taxid <- kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]
adata@meta.data <- cbind(adata@meta.data,kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]) #AddMetaData behaves weirdly!!

#Note, some species NA. example is taxid=1386. listed as Bacillus (Bacillus rRNA group 1), genus, firmicutes

#### How many species?
KrakenSpeciesDistribution(adata)  ### Could use metadata column! TODO   rewrite




################################################################################
################## The rest of KRAKEN analysis #################################
################################################################################


## Dimensional reduction using kraken
DefaultAssay(adata) <- "kraken"

adata <- adata[,adata$nCount_RNA>100] ## Reduce to sensible number
adata

adata <- adata[,adata$species!="Homo sapiens"]  #clearly background!

if(TRUE){
  #ATAC-seq style. think not the best way here
  adata <- RunTFIDF(adata)
  adata <- FindTopFeatures(adata, min.cutoff = 'q0')
  adata <- RunSVD(adata)
  DepthCor(adata)
  adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:30, reduction.name = "kraken_umap")  ## depth seems to be less of a problem here
} else {
  
  adata <- NormalizeData(adata)
  adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
  adata <- ScaleData(adata, features = rownames(adata))
  adata <- RunPCA(adata, features = VariableFeatures(object = adata))
  adata <- RunUMAP(adata, dims = 1:20, reduction.name = "kraken_umap")
}
#DimPlot(object = adata, label = TRUE) + NoLegend()

adata@meta.data$phylum

DimPlot(object = adata, label = TRUE, group.by = "genus", reduction = "kraken_umap")
ggsave("/husky/henriksson/atrandi/wgs_saliva1/umap_genus.pdf", width=15, height=5)

DimPlot(object = adata, label = TRUE, group.by = "species", reduction = "kraken_umap")
ggsave("/husky/henriksson/atrandi/wgs_saliva1/umap_species.pdf", width=40, height=10, limitsize=FALSE)

df <- adata@meta.data
df <- sqldf("select species, count(*) as cnt from df group by species")
df <- df[order(df$cnt),]
df$species <- factor(df$species, levels = df$species)
ggplot(df[df$cnt>5,], aes(species, cnt)) + geom_bar(stat="identity") + coord_flip() + theme_bw()
ggsave("/husky/henriksson/atrandi/wgs_saliva1/hist_species.pdf", width=7, height=10, limitsize=FALSE)




saveRDS(adata, "/husky/henriksson/atrandi/wgs_saliva1/kraken.RDS")


#Compare with depth. "human" cells got few counts and end up in the middle
adata$log_cnt <- log10(1+adata$nCount_RNA)


FeaturePlot(adata, features = "log_cnt")


FeaturePlot(adata, features = "Xanthomonas campestris")

adata$perc_human <- adata@assays$RNA$counts["Homo sapiens",]/colSums(adata@assays$RNA$counts)
FeaturePlot(adata, features = "perc_human")

FeaturePlot(adata, features = "Homo sapiens")


### coverage per species
ggplot(adata@meta.data, aes(species)) + geom_histogram(stat="count") + coord_flip()
df <- adata@meta.data
df <- sqldf("select species, count(*) as cnt, avg(log_cnt) as log_cnt from df group by species")
df <- df[order(df$cnt),]
df



### Compare species kneeplots
KneeplotPerSpecies(adata, max_species = 30)
ggsave("/husky/henriksson/atrandi/wgs_saliva1/kneeplot_species30.pdf", width=7, height=10, limitsize=FALSE)
KneeplotPerSpecies(adata, max_species = 10)
ggsave("/husky/henriksson/atrandi/wgs_saliva1/kneeplot_species10.pdf", width=7, height=10, limitsize=FALSE)



################################################################################
################## Filter out human reads ######################################
################################################################################


### Get reads in fastq format for BWA  --- will not be needed later as the conversion will be made on the fly
BascetMapTransform(
  bascetRoot, 
  "filtered", 
  "asfq",
  out_format="R1.fq.gz"
  #  runner=inst
)

### Perform alignment -- internally wraps mapshard
BascetAlignToReference(
  bascetRoot,
  useReference="/data/henlab/ref_genome/bwa_human/all.fa",
  numLocalThreads=10
)

BascetFilterAlignment(
  bascetRoot,
  inputName="asfq",
  numLocalThreads=10
  
)








################################################################################
################## QUAST todo ##################################################
################################################################################

###
BascetMapCell(
  bascetRoot,
  withfunction = "_quast",
  inputName = "filtered",
  outputName = "quast"
  #  runner=bascet_inst
)

#quast is in /home/mahogny/.local/bin/quast.py

################################################################################
################## Preprocessing for de novo KMERs #############################
################################################################################


### TODO note: now we likely wish to reduce the number of droplets to only cover relevant cells!


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
  runner=bascet_inst
)



### Generate min-hash from each per-cell KMER database
BascetMapCell(
  bascetRoot,
  withfunction = "_minhash",
  inputName = "kmc",
  outputName = "minhash",
  runner=bascet_inst
)





if(TRUE){
  ############################ Could be good to use a #read cutoff; this avoids pulling in noisy kmers. how about weight 1/#read count?
  bascetRoot <- "/husky/henriksson/atrandi/v2_wgs_miseq2"
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
    runner=bascet_inst
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
  runner=bascet_inst)

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
  counts = CreateAssayObject(cnt),
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
################## Reference-based mapping to get "ground truth" ###############
################################################################################


### Perform alignment -- this is a wrapper for mapshard
BascetAlignToReference(
  bascetRoot,
  useReference="/husky/fromsequencer/241016_novaseq_wgs2/trimmed/ref10/all.fa",
  numLocalThreads=10
)

### Generate fragments BED file suited for quantifying reads/chromosome using Signac later -- this is a wrapper for mapshard
BascetBam2Fragments(
  bascetRoot,
  runner=bascet_inst
)


################################################################################
################## Get alignment stats per cell ################################   
################################################################################


#200% cpu. to get more speed, we need to divide up parsing between threads.
#then we need to support addition of matrices
BascetCountChrom(
  bascetRoot, 
  runner=bascet_inst
)


##################### TODO fucked up fragments file!!

### TODO call bascet countchrom


cnt <- ReadBascetCountMatrix("/husky/henriksson/atrandi/wgs_novaseq3/chromcount.1.hd5")  #cnt_al todo delete
#colnames(cnt)
#rownames(cnt)

## TODO: we clearly need 

#adata[["chrom_cnt"]] <- CreateAssayObject(t(cnt)) #TODO swap
adata <- CreateSeuratObject(CreateAssayObject(t(cnt)), assay="chrom_cnt")  #### TODO swap order!!


## Add the counts to our Seurat object
#adata[["chrom_cnt"]] <- FragmentCountsPerChromAssay(bascetRoot)   #### ... first time it crashed on line 300. now fine???
#adata[["chrom_cnt"]] <- FragmentCountsPerChromAssay(bascetRoot, inputName = "new_fragments.1.tsv.gz")   #### ... ignore
DefaultAssay(adata) <- "chrom_cnt"


#Figure out which chromosome has most reads in which cell
cnt <- adata@assays$chrom_cnt$counts
adata$chr_aln <- rownames(cnt)[apply(cnt, 2, which.max)]


adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff = 'q0')
adata <- RunSVD(adata,n = nrow(adata@assays$chrom_cnt$counts)-1) ## need to do this when few chroms
DepthCor(adata)
adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:(nrow(adata@assays$chrom_cnt$counts)-1), reduction.name = "aln_chrom_umap")  ## even 2nd factor is correlated

DimPlot(adata, group.by = "chr_aln", reduction = "aln_chrom_umap", label = TRUE)


# #NC_017316.1   Enterococcus faecalis OG1RF, complete sequence
# #CP086328.1    Bacillus pacificus strain anQ-h4 chromosome, complete genome

#3     CP086328.1 5252926           Bacillus pacificus (ATCC 10987)   big one
#4     CP086329.1  341714           Bacillus pacificus (ATCC 10987)

### why do other species die out if we sum up???

#Log depth is easier to see
adata$log_aln_cnt <- log10(1+colSums(adata@assays$species_cnt$counts))

#Compute purity of genomes
adata$purity_aln <- colMaxs(adata@assays$species_cnt$counts)/colSums(adata@assays$species_cnt$counts)


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
  
  map_seq2strain$id <- stringr::str_replace_all(map_seq2strain$id,stringr::fixed("_"), "-") #### for compatibility with signac
  
  map_seq2strain
}




map_seq2strain <- get_map_seq2strain()[,c("id","len","strain")]
strain_genomesize <- sqldf::sqldf("select sum(len) as len, strain from map_seq2strain group by strain")


########## Produce a count matrix on strain level
DefaultAssay(adata) <- "chrom_cnt"
adata[["species_cnt"]] <- ChromToSpeciesCount(adata, map_seq2strain)  #gives warning. coerce ourselves to dgCMatrix

#Figure out which species has most reads in which cell
cnt <- adata@assays$species_cnt$counts
adata$species_aln <- rownames(cnt)[apply(cnt, 2, which.max)]


#####
adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff = 'q0')
adata <- RunSVD(adata,n = nrow(adata@assays$species_cnt$counts)-1) ## need to do this when few chroms
adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:(nrow(adata@assays$species_cnt$counts)-1), reduction.name = "aln_species_umap")  

DimPlot(adata, group.by = "species_aln", reduction = "aln_species_umap", label = TRUE)

saveRDS(adata, "/data/henlab/temp/aln.RDS")

egg::ggarrange(
  DimPlot(adata, group.by = "species_aln", reduction = "aln_species_umap", label = TRUE),
  DimPlot(adata[,adata$log_aln_cnt>4 & adata$purity_aln>0.9], group.by = "species_aln", reduction = "aln_species_umap"),
  FeaturePlot(adata, features = "log_aln_cnt", reduction = "aln_species_umap"),
  FeaturePlot(adata, features = "purity_aln", reduction = "aln_species_umap")
)
#New species distribution
ggplot(
  adata@meta.data[adata$log_aln_cnt>1 & adata$purity_aln>0.95,],
  aes(species_aln)) + geom_bar() + coord_flip()


KrakenSpeciesDistribution(adata, use_assay = "species_cnt")
### TODO: 
DefaultAssay(adata) <- "species_cnt"
KneeplotPerSpecies(adata)


adata_clean <- adata[,adata$log_aln_cnt>4 & adata$purity_aln>0.9]
adata_clean <- RunTFIDF(adata_clean)
adata_clean <- FindTopFeatures(adata_clean, min.cutoff = 'q0')
adata_clean <- RunSVD(adata_clean,n = nrow(adata_clean@assays$species_cnt$counts)-1) ## need to do this when few chroms
adata_clean <- RunUMAP(object = adata_clean, reduction = 'lsi', dims = 1:(nrow(adata_clean@assays$species_cnt$counts)-1), reduction.name = "aln_species_umap")  


DimPlot(adata_clean, group.by = "species_aln", reduction = "aln_species_umap")
ggsave("/home/mahogny/species_aln.pdf", width = 10, height = 6) #### for february presentation




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
BarnyardPlotMatrix(adata[
  1:2,
  adata$log_aln_cnt>4])


subdata <- adata[
  c("Streptococcus mutans (ATCC 700610)","Clostridium beijerinckii (ATCC 35702)"),
  adata$log_aln_cnt>4]

df <- data.frame(
  x=adata@assays$species_cnt$counts[1,],
  y=adata@assays$species_cnt$counts[2,]
)
ggplot(df, aes(x,y)) + geom_point() + scale_x_log10() + scale_y_log10() 

#
BarnyardPlotMatrix()


rownames(adata)

## TODO: after detecting doublets, we should color by this in the plot   ; color by any metadata!!





################################################################################
################## Correlation of genomes ######################################
################################################################################

DefaultAssay(adata) <- "species_cnt"


SpeciesCorrMatrix(adata)


# other code in /home/mahogny/jupyter/scWGS









get_idx_stats_from_bam <- function(fname){
  idx <- read.table(pipe(paste("samtools idxstats ",fname)),sep="\t")
  colnames(idx) <- c("id","len","mapped","unmapped")
  idx  
}
idxstat <- get_idx_stats_from_bam("/husky/henriksson/atrandi/wgs_saliva1/sorted_al/aligned.1.bam")


### 98% unmapped
idxstat$unmapped[idxstat$id=="*"]/(sum(idxstat$unmapped)+sum(idxstat$mapped))








###########################


bascet_runner.default <- LocalRunner(direct = TRUE, show_script=FALSE)
#bascet_instance.default <- BascetInstance(bin = "/home/mahogny/github/bascet/target/debug/bascet", tempdir = "/tmp/")
bascet_instance.default <- getBascetSingularityImage("/home/mahogny/github/bascet/singularity/") #BascetInstance(bin = "/home/mahogny/github/bascet/target/debug/bascet", tempdir = "/tmp/")#, direct = TRUE, show_script=FALSE)

GetDefaultBascetInstance()

bascetRoot <- "/husky/henriksson/atrandi/v2_wgs_miseq2"
bascetRoot <- "/husky/henriksson/atrandi/wgs_novaseq3"
all_kmer <- AggregateMinhashes(bascetRoot) 


streamer <- extractstreamer_start(bascet_instance = bascet_instance)
extractstreamer_open(streamer, "/husky/henriksson/atrandi/v2_wgs_miseq2/minhash.1.zip")


