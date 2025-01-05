source("somezorn.R")

################################################################################
################ Zorn API for running Bascet via SLURM #########################
################################################################################


# These commands take care of producing job arrays

# also check https://github.com/USCbiostats/slurmR


# BascetSettings(path="~/bascetbin")  #can be passed if needed
# LocalInstance()  #can be given instead, if running locally


####### Where all bascets are stored
bascetRoot <- "/husky/cellbuster/randomproj"


####### Default SLURM instance settings
slurmsettings <- SlurmInstance(partition="shared", account="blet", time="0-24:00:00")


####### Generate BAM with barcodes from input raw FASTQ
rawmeta <- DetectRawFileMeta("/husky/fromsequencer/241206_novaseq_wgs3/raw")
p10 <- BascetGetRawAtrandiWGS(
  bascetRoot, 
  rawmeta, 
  runner=SlurmInstance(slurmsettings, ncpu=5, mem="5g")
)
WaitForJob(p10)  #Has possibility of ctrl+c; just keeps polling, possibly with a status indicator from log. or keep plotting log file

#possibly useful:
# CancelJob(p1)  
# CancelAllSlurmJobs()
# JobStatus(p1)
# JobLog(p1)


####### Generate FASTQ ZIP file
p0 <- BascetPartition(
  bascetRoot,
  runner=SlurmInstance(slurmsettings, ncpu=10, mem="60g")
)
WaitForJob(p0) 


####### Generate zip file with fastq from each cell
# useful if this command should be able to take multiple debarcoded files as input; for now, merge with samtools before
p1 <- BascetPartition(
  bascetRoot,
  ## default settings assumed here
  runner=SlurmInstance(slurmsettings, ncpu=10, mem="60g")
)
WaitForJob(p1)  



####### Run Spades to assemble the genomes  -- move to SKESA
p2 <- BascetAssemble(
  bascetRoot,
  runner=SlurmInstance(slurmsettings, ncpu=10, mem="60g")
)
WaitForJob(p2)


####### Build kmer database
p3 <- BascetCount(
  bascetRoot,
  runner=SlurmInstance(slurmsettings, ncpu=10, mem="10g")
)
WaitForJob(p3)



####### Select kmers that appear useful for clustering
p4 <- BascetFeaturise( 
  bascetRoot,
  runner=SlurmInstance(slurmsettings, ncpu=10, mem="10g")
)
WaitForJob(p4)



####### Build count table from kmer table and selected kmers
p5 <- BascetQuery(
  bascetRoot,
  runner=SlurmInstance(slurmsettings, ncpu=5, mem="5g")
)
WaitForJob(p5)




################################################################################
################ Extra things to support later #################################
################################################################################


####### Add raw FASTQs of isolates, enabling them to be treated as cells, assembled etc
# todo think about what to name this file. maybe separate from cell fastq to enable easy rerunning
# todo allow this function to be called multiple times?
pxx <- BascetAddIsolateRawFastq(
  bascetRoot,
  listFastqR1 = c(...),
  listFastqR2 = c(...),
  names = c(...), ### give each isolate a more sane name. default is name of fastq otherwise 
  runner=SlurmInstance(slurmsettings, ncpu=5, mem="5g")
)



####### Add isolate genomes, treating them as assembled, enabling clustering, comparison with cells, etc
# todo think about what to name this file. maybe separate from cell assembly to enable easy rerunning
# todo allow this function to be called multiple times?
pxx <- BascetAddAssembledIsolate(
  bascetRoot,
  listFasta = c(...),
  names = c(...), ### give each isolate a more sane name. default is name of fasta otherwise 
  runner=SlurmInstance(slurmsettings, ncpu=5, mem="5g")
)


################################################################################
################ Additional QC #################################################
################################################################################



######### Call Quast for all cells
p6 <- BascetMapCell(
  bascetRoot,
  "quast",
  runner=SlurmInstance(slurmsettings, ncpu=5, mem="5g")
)

WaitForJob(p6)




######### Example callback function for aggregating data
aggr.quast <- function(bascetFile, cellID){
  ###### Option #1
  tmp <- BascetReadMapFile(bascetFile, cellID, "out.csv", as="tempfile")
  my_df <- read.csv(tmp)
  file.remove(tmp)

  ###### Option #2
  my_data <- BascetReadMapFile(bascetFile, cellID, "out.csv", as="text")   #equivalent to readLines() i.e. one big string returned (or get a list, one per line?)

  
  return(data.frame(
    quality=666,
    completeness=50
  ))
}





######### Aggregate data from previous Map call
quast_aggr <- MapListAsDataFrame(BascetAggregateMap(
  bascetRoot,
  "quast",
  aggr.quast
  ### Option: can later add use.runner=SlurmInstance(slurmsettings, ncpu=10, mem="60g")   to run this in parallel via slurm
))







################# robert commands beneath
#robert partition -i file_in
#robert assemble
#robert count
#robert featurise
#robert query



################################################################################
################ Loading count table, processing etc ###########################
################################################################################





######### Example callback function for aggregating data
aggr.quast <- function(bascetFile, cellID){
  ###### Option #1
  tmp <- BascetReadMapFile(bascetFile, cellID, "out.csv", as="tempfile")
  my_df <- read.csv(tmp)
  file.remove(tmp)

  ###### Option #2
  my_data <- BascetReadMapFile(bascetFile, cellID, "out.csv", as="text")   #equivalent to readLines() i.e. one big string returned (or get a list, one per line?)

  
  return(data.frame(
    quality=666,
    completeness=50
  ))
}





LoadBascetCounts <- function(bascetRoot, bascetFile="somedefault"){
  
}





cnt <- LoadBascetCounts(bascetRoot)


#New Seurat object etc










################################################################################
################ latest testing ################################################
################################################################################


## Cells having contigs
BascetCellNames("/home/mahogny/jupyter/bascet/zorn","contigs")





