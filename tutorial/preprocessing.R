library(zorn)
library(ggplot2)


### Decide how to run Bascet. Currently it is set up to run locally, synchronized (easiest for beginners)
inst <- LocalInstance(direct = TRUE, showScript=TRUE)


### Decide where to store all processed data. This directory must exist!
bascetRoot = "testdata"  


################################################################################
################## Sort reads per cell etc #####################################
################################################################################

### Figure out what input read files to use. This method detects this automatically
### based on file names
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

### Shardify selected cells i.e. divide into multiple sets of files for parallel processing
BascetShardify(
  bascetRoot,
  includeCells = includeCells,
  runner = inst
)




################################################################################
################## One approach: KRAKEN preproceesing ##########################
################################################################################



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
################## One approach: Min-hash preproceesing ########################
################################################################################

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



### Build count table by looking up selected KMERs in per-cell KMER databases
BascetQuery(
    bascetRoot, 
    useKMERs = useKMERs,
    runner=inst)








