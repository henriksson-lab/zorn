library(zorn)

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

CreateKmerAssay <- function(counts) {
  chrom_assay <- CreateAssayObject(  ##################### in the future, can add other metadata in here too for visualization?
    counts = counts
  )
  chrom_assay
}


adata <- CreateSeuratObject(
  counts = CreateKmerAssay(cnt),
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









### Get reads in fastq format
BascetMapTransform(
  bascetRoot, 
  "filtered", 
  "for_bwa",
  out_format="fq.gz",   ## but we need two fq as out!! ideally at least. or if R1.fq.gz => write two of them. otherwise gather?
  runner=inst)



### Perform alignment  
#.... where to keep scripts? this is a pure-zorn feature? unless we support piping formats in and out on the fly (i.e. transform to fastq included)
#.... putting this in bascet only makes sense if bascet can pipe input directly to the receiving software. and out of it! 
BascetMapShard(
  bascetRoot,
  withfunction = "_bwa",
  inputName = "for_bwa",
  outputName = "aligned",
  runner=inst
)







cmd <- paste("bwa mem",
             "/husky/fromsequencer/241016_novaseq_wgs2/trimmed/ref10/all.fa",
             file.path(bascetRoot,"for_bwa.1.fq.gz"),
             "-t 10",
             "> ",
             file.path(bascetRoot,"aligned.1.bam")
)
system(cmd)


## Check output
system(
  paste("samtools view", file.path(bascetRoot,"aligned.1.bam"))
)

## Produce bed-file
system(
  paste("samtools view", file.path(bascetRoot,"aligned.1.bam"))
)


system(
  paste(
    "bedtools bamtobed -i ",
    file.path(bascetRoot,"aligned.1.bam"),
    " > ",
    file.path(bascetRoot,"aligned_bed.1.bam")
  )
)
# NZ_CP009792.1	1245738	1245889	BASCET_F3_H5_B8_D11::71	60	-

bedtable <- read.table(file.path(bascetRoot,"aligned_bed.1.bam"),sep="\t")
colnames(bedtable) <- c("chr","from","to","cell_id","score","strand")
bedtable$cell_id <- stringr::str_remove(bedtable$cell_id,"BASCET_")
bedtable$cell_id <- stringr::str_split_i(bedtable$cell_id,":",1)
bedtable <- bedtable[order(bedtable$chr, bedtable$from, bedtable$to),]



j### Format as 10x BED file
#chr1	10067	10333	CGATCCTTCGCTCCAT-1	1
as_fragments <- data.frame(
  chr=bedtable$chr,
  from=bedtable$from,
  to=bedtable$to,
  cell_id=bedtable$cell_id,
  stupid=1
)
fname_fragments <- file.path(bascetRoot,"fragments.tsv")
write.table(as_fragments, fname_fragments, sep="\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

system(paste("bgzip -f",fname_fragments))
system(paste0("tabix -p bed ",fname_fragments,".gz"))


################################################################################
################## Get alignment stats per cell ################################
################################################################################

#Table of cell_id vs chromosome count
chr_cnt <- sqldf::sqldf("select chr, cell_id, count(*) as cnt from bedtable group by chr, cell_id")
chr_cnt <- reshape2::acast(chr_cnt, cell_id~chr, value.var = "cnt", fill = 0)

chr_cnt

#Align to adata order
chr_cnt <- chr_cnt[colnames(adata),]
which.max(chr_cnt)

adata$max_chr <- colnames(chr_cnt)[apply(chr_cnt, 1, which.max)]

DimPlot(adata, group.by = "max_chr")

#NC_017316.1   Enterococcus faecalis OG1RF, complete sequence
#CP086328.1    Bacillus pacificus strain anQ-h4 chromosome, complete genome



################################################################################
################## Use Signac as feature counter ###############################
################################################################################






