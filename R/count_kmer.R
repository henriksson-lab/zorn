###############################################
############################################### de novo kmer system
###############################################



###############################################
#' Compute minhashes for each cell.
#' This is a thin wrapper around BascetMapCell
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param maxReads The maximum number of reads per cell to sample
#' @param kmerSize The KMER size for the hashing
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetComputeMinhash <- function( 
    bascetRoot, 
    inputName="filtered", 
    outputName="minhash", 
    overwrite=FALSE,
    maxReads=100000,  #This is most likely enough to get an overall histogram
    kmerSize=31,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  #check arguments
  stopifnot(is.integer.like(kmerSize))
  stopifnot(is.integer.like(maxReads))
  
  #most arguments checked in thgis call
  BascetMapCell(
    bascetRoot=bascetRoot, 
    withfunction="_minhash_fq", 
    inputName=inputName, 
    outputName=outputName,
    #includeCells=includeCells
    args = list(
      KMER_SIZE=format(kmerSize, scientific=FALSE),
      MAX_READS=format(maxReads, scientific=FALSE)
    ),
    overwrite=overwrite,
    runner=runner,
    bascetInstance=bascetInstance)
}




###############################################
#' Gather all minhashes into a single histogram file
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard (should contain minhashes)
#' @param outputName Name of output file
#' @param includeCells List of cells to include, or NULL for all cells
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetMakeMinhashHistogram <- function( 
  bascetRoot, 
  inputName="minhash", 
  outputName="minhash_hist.csv", 
  includeCells=NULL,
  overwrite=FALSE,
  runner=GetDefaultBascetRunner(),
  bascetInstance=GetDefaultBascetInstance()
){
  #check arguments
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.valid.listcells(includeCells))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance) 
  
  #Figure out input and output file names
  inputFiles <- file.path(bascetRoot, detectShardsForFile(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }
  
  outputFile <- file.path(bascetRoot, outputName)

  #If cell list is provided, produce a file for input (not all transform calls can handle this, so optional)
  produce_cell_list <- !is.null(includeCells)
  if(produce_cell_list) {
    #Currently using the same cell list for all shards (good idea?)
    list_cell_for_shard <- list()
    for(i in 1:length(inputFiles)){
      list_cell_for_shard[[i]] <- includeCells
    }
  }

  if(bascetCheckOverwriteOutput(outputFile, overwrite)) {
    #Make the command
    cmd <- c(
      #shellscript_set_tempdir(bascetInstance),
      if(produce_cell_list) shellscriptMakeFilesExpander("CELLFILE", list_cell_for_shard),
      paste(
        bascetInstance@prependCmd,
        bascetInstance@bin, 
        "minhash-hist",
        if(produce_cell_list) "--cells ${CELLFILE[$TASK_ID]}",
        "-t $BASCET_TEMPDIR",
        "-i", shellscriptMakeCommalist(inputFiles),
        "-o", outputFile
      )
    )
  
    #Run the job
    RunJob(
      runner = runner, 
      jobname = "bascet_minhash_hist",
      bascetInstance = bascetInstance,
      cmd = cmd,
      arraysize = 1
    )  
  } else {
    new_no_job()
  }
}








###############################################
#' Build count table from FASTQ reads and a list of selected kmers
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param useKMERs List of KMERs to query
#' @param maxReads The maximum number of reads per cell to sample
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetQueryFq <- function( #666
    bascetRoot, 
    inputName="filtered",
    outputName="kmer_counts", 
    useKMERs,
    maxReads=1000000, 
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  #check arguments
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  #useKMERs TODO
  stopifnot(is.positive.integer(maxReads))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance) 
  
  
  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, detectShardsForFile(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }
  
  if(length(useKMERs)==0){
    stop("No KMERs")
  }
  
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "h5", num_shards)
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    #Run the job
    RunJob(
      runner = runner, 
      jobname = "bascet_query",
      bascetInstance = bascetInstance,
      cmd = c(
        #shellscript_set_tempdir(bascetInstance),
        shellscriptMakeBashArray("files_in",inputFiles),
        shellscriptMakeBashArray("files_out",outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),
        
        shellscriptMakeOneFileExpander("KMERFILE", useKMERs), 
        paste(
          bascetInstance@prependCmd,
          bascetInstance@bin, 
          "query-fq",
          "-t $BASCET_TEMPDIR",
          "-m ", format(maxReads, scientific=FALSE),
          "-f $KMERFILE",
          "-i ${files_in[$TASK_ID]}",
          "-o ${files_out[$TASK_ID]}"
        )
      ),
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
}




###############################################
#' Read histogram of KMERs, the output of BascetMakeMinhashHistogram
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' 
#' @return KMER histogram
#' @export
BascetReadMinhashHistogram <- function(
    bascetRoot,
    inputName="minhash_hist.csv"
){
  #check arguments
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.character(inputName))

  #Read file
  fname <- file.path(bascetRoot,inputName)
  if(!file.exists(fname)) {
    stop(paste("File is missing:",fname))
  }
  dat <- as.data.frame(data.table::fread(fname))
  colnames(dat) <- c("kmer","cnt")
  dat <- dat[order(dat$cnt, decreasing = TRUE),,drop=FALSE]
  dat
}









###############################################
#' Pick random KMERs from KMC3 database. The choice is among KMERs within a frequency range
#' 
#' @param kmerHist Data.frame with KMER frequency
#' @param minFreq Pick KMERs having minimum frequency
#' @param maxFreq Pick KMERs having maxiumum frequency
#' 
#' @return List of KMERs
#' @export
ChooseInformativeKMERs <- function(
    kmerHist,
    #numPick=1000, 
    minFreq=0.005, 
    maxFreq=1
) {
  #check arguments
  stopifnot(is.data.frame(kmerHist))
  stopifnot(is.numeric.range01(minFreq))
  stopifnot(is.numeric.range01(maxFreq))
  stopifnot(minFreq < maxFreq)
  
  #maxFreq = 1
  #minFreq = 0.02
  
  ## Possibly expensive to get them all... is the total count stored somewhere? can we sample?
  kmerHist <- kmerHist[order(kmerHist$cnt),]
  num_kmer <- nrow(kmerHist)
  
  max_cnt <- max(kmerHist$cnt)
  
  #all_kmer <- kmerHist$kmer[round(num_kmer*minFreq):round(num_kmer*maxFreq)]
  sub_kmerHist <- kmerHist[kmerHist$cnt>=max_cnt*minFreq & kmerHist$cnt<=max_cnt*maxFreq,]
  #[kmerHist$cnt>max_cnt*minFreq & kmerHist$cnt<max_cnt*maxFreq]
  
  #count_range <- kmerHist$cnt[c(round(num_kmer*minFreq),round(num_kmer*maxFreq))]
  print(paste(
    "Counts will be in range ",
    min(sub_kmerHist$cnt), 
    max(sub_kmerHist$cnt),
    "#kmers", nrow(sub_kmerHist)))

  all_kmer <- sub_kmerHist$kmer
  
  #if(length(all_kmer)<numPick){
  #  stop(paste("Cannot pick",numPick," KMERS; only",length(all_kmer),"pass frequency criteria"))  
  #}
  
  all_kmer
  #set.seed(0)
  #sample(all_kmer, numPick)
}

