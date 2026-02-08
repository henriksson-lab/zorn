

################################################################################
################ Bascet command line tool: debarcoding raw fastq ###############
################################################################################



###############################################
#' Detect metadata for raw input FASTQ files
#' 
#' _R1 -- common from illumina sequencer
#' SRR****_1.fastq.gz -- typical from SRA
#' 
#' TODO Would be convenient to handle multiple samples, as sample1/xxx; in this case, 
#' should prepend the sample name to the barcodes.
#' 
#' issue: when shardifying, good to keep info about what to merge. this reduces the work plenty!
#' could keep a list of which shards belong for the next step
#' 
#' 
#' @param rawRoot Path to folder with FASTQ files
#' @param verbose Print additional information, primarily to help troubleshooting
#' 
#' @return A data frame with metadata for the raw input files
#' @export
DetectRawFileMeta <- function(
    rawRoot, 
    verbose=FALSE
){
  #check arguments
  stopifnot(dir.exists(rawRoot))
  stopifnot(is.logical(verbose))
  
  
  #rawRoot <- "/husky/fromsequencer/241206_novaseq_wgs3/raw"
  allfiles <- list.files(rawRoot)
  
  if(verbose){
    print("Found the following files in the directory")
    print(allfiles)
  }
  
  allfiles <- allfiles[
    stringr::str_ends(allfiles,".fq.gz") | 
      stringr::str_ends(allfiles, ".fastq.gz") | 
      stringr::str_ends(allfiles, ".fq") | 
      stringr::str_ends(allfiles,".fastq")]
  
  if(verbose){
    print("Found the following raw read-like files")
    print(allfiles)
  }
  
  r1_files <- allfiles[stringr::str_detect(allfiles,stringr::fixed("_R1"))]
  
  if(verbose){
    print("Found the following R1-like files")
    print(r1_files)
  }
  
  ## TODO for SRA, replace with _1.fastq.gz
  r2_corresponding <- stringr::str_replace(r1_files,stringr::fixed("R1"),"R2")
  
  if(!all(r2_corresponding %in% allfiles)){
    stop("Not all R1 files appear to have a corresponding R2 file")
  }
  
  ### Gather first set of metadata
  meta <- data.frame(
    r1=r1_files,
    r2=r2_corresponding,
    dir=file.path(rawRoot) #this standarizes the trailing /
  )
  meta$prefix <- ""
  
  #Guess prefixes
  meta$possible_prefix <- stringr::str_split_i(meta$r1, "_S[0123456789]",1)
  unique_prefix <- unique(meta$possible_prefix)
  if(verbose){
    print("Detected possible prefixes:")
    print(unique_prefix)
  }
  
  #If there is more than one prefix, then we have to add them. otherwise just keep it simple
  if(length(unique_prefix)>1){
    #Sanitize prefixes. Some characters will break BAM tags etc
    meta$prefix <- meta$possible_prefix
    meta$prefix <- stringr::str_remove_all(meta$prefix, " ")
    meta$prefix <- stringr::str_remove_all(meta$prefix, "/")
    meta$prefix <- stringr::str_remove_all(meta$prefix, "\"")
    
    print("Detected multiple libraries")    
  } else {
    print("Detected a single library")    
    
    if(verbose){
      print("Detected only one possible prefix, so not adding it")
    }
  }
  
  #Return metadata
  meta[,c("prefix","r1","r2","dir")]
}

#rawmeta <- DetectRawFileMeta("/husky/fromsequencer/241206_novaseq_wgs3/raw")




###############################################
#' Extract barcodes and trim input raw FASTQ
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param rawmeta Metadata for the raw FASTQ input files. See DetectRawFileMeta
#' @param maxShardSize Estimated maximum size of output shard. Can be set higher but as sorting is also performed during sharding, it can be overall more efficient to only do partial sorting during this command
#' @param outputName Name output files: Debarcoded reads
#' @param chemistry The type of data to be parsed
#' @param barcodeTolerance Optional: Number of mismatches allowed in the barcode for it to still be considered valid
#' 
#' @param numThreads Number of threads to use per job. Default is the number from the runner
#' @param numReadThreads Advanced setting: Number of threads for reading 
#' @param numDebarcodeThreads Advanced setting: Number of threads for debarcoding
#' @param numSortingThreads Advanced setting: Number of threads for sorting, first phase
#' @param numMergeSortingThreads Advanced setting: Number of threads for sorting, second phase
#' @param numWriteThreads Advanced setting: Number of threads for writing
#' @param numCompressThreads Advanced setting: Number of threads for compressing
#' 
#' @param totalMem Total memory to allocate
#' @param streamBufferSize Advanced setting: Stream buffer size (fraction, given as e.g. "10%")
#' @param sortBufferSize Advanced setting: Sort buffer size (fraction, given as e.g. "10%")
#' @param pageBufferSize Advanced setting: Page buffer size (fraction, given as e.g. "10%")
#' @param compressionLevel Advanced setting: Compression level (0..12)
#' 
#' @param overwrite 
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetGetRaw <- function(
    bascetRoot, 
    rawmeta,
    maxShardSize="200g",  ### if any?   maybe no need!
    outputName="debarcoded", 
    chemistry=c("atrandi-wgs","atrandi-wgslr","atrandi-rnaseq","parse-bio"),  #TODO any way to get list from software?
    subchemistry=NULL,
    barcodeTolerance=NULL,
    
    numThreads=NULL,
    numReadThreads=NULL,
    numDebarcodeThreads=NULL,
    numSortingThreads=NULL,
    numMergeSortingThreads=NULL,
    numWriteThreads=NULL,
    numCompressThreads=NULL,
    
    totalMem=NULL,
    streamBufferSize=NULL,
    sortBufferSize=NULL,
    compressBufferSize=NULL,
    compressRawBufferSize=NULL,
    
    compressionLevel=NULL,

    numMergeStreams=NULL,

    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  #Ensure there are files to process
  if(is.data.frame(rawmeta) && nrow(rawmeta)==0){
    stop("No input files")
  }
  
  #Set number of threads if not given
  if(is.null(numThreads)) {
    numThreads <- as.integer(runner@ncpu)
  }
  
  #Need a minimum number of threads
  if(numThreads < 6) {
    print("Note: Setting number of threads to 6 as this is the minimum, even if the CPU has fewer cores")
    numThreads <- 6
  }
  
  #Check input arguments 
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.data.frame(rawmeta))
  stopifnot(is.valid.shardname(outputName))
  chemistry <- match.arg(chemistry)
  stopifnot(is.numeric(barcodeTolerance) || is.null(barcodeTolerance))
  stopifnot(is.valid.threadcount(numThreads))

  stopifnot(is.null(numMergeStreams) || is.positive.integer(numMergeStreams))
  if(!is.null(numMergeStreams)) {
    stopifnot(is.positive.integer(numMergeStreams))
    numMergeStreams <- as.integer(numMergeStreams)
    stopifnot(numMergeStreams >= 2)
  }

  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))

  #Check compression setting  
  if(!is.null(compressionLevel)) {
    stopifnot(is.integer(compressionLevel))
    compressionLevel <- as.integer(compressionLevel)
    stopifnot(compressionLevel >= 0 && compressionLevel <= 12)
  }
  
  #Check memory sizes
  stopifnot(is.null(streamBufferSize) || is.percent.string(streamBufferSize))
  stopifnot(is.null(sortBufferSize) || is.percent.string(sortBufferSize))
  stopifnot(is.null(compressBufferSize) || is.percent.string(compressBufferSize))
  stopifnot(is.null(compressRawBufferSize) || is.percent.string(compressRawBufferSize))
  if(!is.null(totalMem)) {
    totalMem <- parse_size_string(totalMem)
    stopifnot(totalMem > fs::fs_bytes("1Gb"))
  } else {
    #Take memory from runner if possible
    if(runner@mem!="") {
      totalMem <- parse_size_string(runner@mem) - fs::fs_bytes(bascetInstance@containerMem)
      stopifnot(totalMem > fs::fs_bytes("1Gb"))
    } else {
      print("Warning: Total memory was not specified. We strongly encourage doing this to ensure performance")
    }
  }

  #Convert size to bytes and check argument
  maxShardSize <- parse_size_string(maxShardSize)

  #Figure out how many output files are needed.
  #Do this by checking size of input files
  rawmeta$filesize <- NA
  for(i in 1:nrow(rawmeta)){
    size_r1 <- file.info(file.path(rawmeta$dir, rawmeta$r1[i]))$size
    size_r2 <- file.info(file.path(rawmeta$dir, rawmeta$r2[i]))$size
    rawmeta$filesize[i] <- fs::as_fs_bytes(size_r1+size_r2)
  }
  rawmeta$need_num_outputs <- ceiling(rawmeta$filesize/maxShardSize)
  

#print(6666)
#print(rawmeta)
#print(rawmeta$need_num_outputs)

  if(any(is.na(rawmeta$filesize))){
    print(rawmeta)
    stop("Not all input files exist or all accessible")
  }
  #  print(rawmeta)
  
  
  # Figure out output file names
  num_shards <- sum(rawmeta$need_num_outputs)
  
  #  print(rawmeta)
  #  print(num_shards)
  
  outputFilesComplete <- makeOutputShardNames(bascetRoot, outputName, "tirp.gz", num_shards)
  
  #Check if libnames should be added
  add_libnames <- any(rawmeta$prefix!="")
  
#print("---")
#print("---")

  #Duplicate rawmeta for each output file
  rawmeta_tostore <- NULL
  cur_start_shard <- 0
  arg_outputFiles <- NULL
  for(i in 1:nrow(rawmeta)) {
    rawmeta_one <- rawmeta[i,,drop=FALSE]
    cur_ids <- cur_start_shard + (1:rawmeta_one$need_num_outputs)
#    print(cur_ids)
    rawmeta_one <- rawmeta_one[rep(1, length(cur_ids)),,drop=FALSE] #replicate rows
    rawmeta_one$shard <- cur_ids
    rownames(rawmeta_one) <- NULL
#print(rawmeta_one)
    rawmeta_tostore <- rbind(rawmeta_tostore, rawmeta_one)
    
    #Also create a ,-separated file list for running Bascet
    arg_outputFiles <- cbind(
      arg_outputFiles,
      stringr::str_flatten(outputFilesComplete[cur_ids], collapse = ",")
    )
    
    #Move to next set of outputs
    cur_start_shard <- cur_start_shard + rawmeta_one$need_num_outputs
  }

  print(rawmeta_tostore)
  #  print(rawmeta)
  #  print(arg_outputFiles)
  #  return(666)
  
  if(bascetCheckOverwriteOutput(outputFilesComplete, overwrite)) {
    
    #Write a file describing the libraries. Store the long format.
    #Important to only write a new file if we also run the jobs.
    #It would be nice to be even more sure that output is properly synched
    write.csv(
      rawmeta_tostore,
      file=file.path(bascetRoot, paste0(outputName, ".meta")),
      row.names = FALSE
    )
    
    RunJob(
      runner = runner, 
      jobname = "Zgetraw",
      bascetInstance = bascetInstance,
      cmd = c(
        #shellscript_set_tempdir(bascetInstance),
        shellscriptMakeBashArray("files_r1",file.path(rawmeta$dir, rawmeta$r1)),
        shellscriptMakeBashArray("files_r2",file.path(rawmeta$dir, rawmeta$r2)),
        shellscriptMakeBashArray("libnames",rawmeta$prefix),
        shellscriptMakeBashArray("files_out",arg_outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),
        
        assembleBascetCommand(bascetInstance, c(
          "get-raw",
          "--temp=$BASCET_TEMPDIR",
          
          paste0("--threads=", numThreads),
          if(!is.null(numReadThreads)) paste0("--countof-threads-read=",numReadThreads),
          if(!is.null(numDebarcodeThreads)) paste0("--countof-threads-debarcode=",numDebarcodeThreads),
          if(!is.null(numSortingThreads)) paste0("--countof-threads-sort=",numSortingThreads),
          if(!is.null(numWriteThreads)) paste0("--countof-threads-write=",numWriteThreads),
          
          if(!is.null(totalMem)) paste0("--memory=",format_size_bascet(totalMem)), 
          if(!is.null(streamBufferSize)) paste0("--sizeof-stream-buffer=",streamBufferSize),
          if(!is.null(sortBufferSize)) paste0("--sizeof-sort-buffer=",sortBufferSize), 
          if(!is.null(compressBufferSize)) paste0("--sizeof-compress-buffer=",compressBufferSize),
          if(!is.null(compressRawBufferSize)) paste0("--sizeof-compress-raw-buffer=",compressRawBufferSize),

          if(!is.null(numMergeStreams)) paste0("--countof-merge-streams=",numMergeStreams),

          if(!is.null(compressionLevel)) paste0("--compression-level=",compressionLevel),
          
          "--r1=${files_r1[$TASK_ID]}",
          "--r2=${files_r2[$TASK_ID]}",
#########          if(add_libnames) "--libname=${libnames[$TASK_ID]}",
          "--out=${files_out[$TASK_ID]}",                 #Each job produces a single output
          chemistry,
          if(!is.null(subchemistry)) paste0("--subchemistry=",subchemistry)
        ))
      ),
      arraysize = nrow(rawmeta)
    )    
  } else {
    new_no_job()
  }
}





###############################################
#' Prepare to shard reads by collecting statistics about
#' each barcode, and filtering out cells with few reads
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param minQuantile Read count-based cutoff for inclusion in final shards
#' @param bascetInstance A Bascet instance
#' @param verbose Print additional information, primarily to help troubleshooting
#' 
#' @return Statistics about the debarcoded reads
#' @export
PrepareSharding <- function(
    bascetRoot, 
    inputName="debarcoded", 
    minQuantile=0.5,
    bascetInstance=GetDefaultBascetInstance(),
    verbose=TRUE
){
  #Check arguments 
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.numeric.range01(minQuantile))
  stopifnot(is.bascet.instance(bascetInstance))
  stopifnot(is.logical(verbose))

  
  #Read metadata to know which files need to be concatenated.
  #This speeds up the process
  fstats <- file.path(bascetRoot, paste0(inputName,".meta"))
  if(!file.exists(fstats)) {
    stop("Missing .meta file for input")
  }
  meta <- read.csv(fstats)
  meta[is.na(meta)] <- ""
  #print(meta)
  meta$prefix <- as.factor(meta$prefix)
  meta$group <- as.integer(meta$prefix)
  group_prefix <- levels(meta$prefix)
  numgroup <- max(as.integer(meta$prefix))
  
  
  #print(meta)
  #print(group_prefix)
  
  #Get all the TIRPs, sum up the reads  
  #OPTION: store shard file name in meta. check that they exist!!
  inputFiles <- detectShardsForFile(bascetRoot, inputName)
  inputFiles <- naturalsort::naturalsort(inputFiles)

  if(length(inputFiles) != nrow(meta)){
    print("=== Detected input files: ===")
    print(inputFiles)
    stop(paste0("Metadata has different number of files (",nrow(meta),") than currently detected (",length(inputFiles),")"))
  }
  
  if(verbose) {
    print("Loading barcode counts")
    pb <- txtProgressBar(min = 0, max = length(inputFiles), initial = 0) 
  }
  
  ### Read BC info for each file
  list_hist <- list()
  for(cur_i in 1:length(inputFiles)) {
    hist_f <- paste0(file.path(bascetRoot, inputFiles[cur_i]),".hist") ## only support TIRP for now
    dat <- readr::read_tsv(  ### use some data.table alternative?
      file = hist_f,
      col_names = FALSE,
      col_types = list(
        #bc = 
          readr::col_character(),
        #cnt = 
          readr::col_double()
      ),
      progress=FALSE
    )
    colnames(dat) <- c("bc","cnt")
    list_hist[[cur_i]] <- dat
    
    if(verbose){
      setTxtProgressBar(pb,cur_i)
    }
  }
  if(verbose) {
    close(pb)
  }  
  
  ### Perform barcode merging
  if(verbose) {
    print("Merging barcode lists for each library")
    pb <- txtProgressBar(min = 0, max = numgroup, initial = 0) 
  }
  list_group_bc_count <- list()
  for(g in 1:numgroup) {

    #Gather counts for cells across files for each library (e.g. multiple lanes)
    dat <- data.table::rbindlist(list_hist[meta$shard[meta$group==g]])
    colnames(dat) <- c("cellid","onecount")
    
    if(verbose){
      #TODO ensure data.table subsetting is ok in library
      print("Example counts")
      print(head(dat))
    }

    #Add up cells across barcodes
    sumdat <- dat[ ,list(count=sum(onecount)), by=cellid]
    
    list_group_bc_count[[g]] <- sumdat[order(count,decreasing=TRUE)]
    
    if(verbose){
      setTxtProgressBar(pb,g) ### or?
    }
  }
  if(verbose) {
    close(pb)
  }  
  
  
  ### Kneeplot analysis for each library
  if(verbose) {
    print("Kneeplot analysis for each library")
    pb <- txtProgressBar(min = 0, max = numgroup, initial = 0) 
  }
  list_knee <- list()
  list_picked_cells <- list()
  for(g in 1:numgroup) {
    #if(verbose) {
    #  print(g)
    #}
    
    #Produce a knee plot. Subsample to speed up the code (this might cause bias?)
    knee <- list_group_bc_count[[g]]
    knee$index <- 1:nrow(knee)
    
    sub_knee <- knee[sample(1:nrow(knee), min(nrow(knee),10000)),, drop=FALSE]
    
    sub_knee$cs <- cumsum(sub_knee$count)
    sub_knee$picked <- sub_knee$cs > sum(sub_knee$count)*minQuantile
    sub_knee$group <- g
    
    #min_quantile <- 0.99
    read_cutoff <- quantile(sub_knee$count, minQuantile)
    sub_knee$picked <- sub_knee$count > read_cutoff
    
    list_knee[[g]] <- sub_knee
    list_picked_cells[[g]] <- knee$cellid[knee$count > read_cutoff]
    
    if(verbose){
      setTxtProgressBar(pb,g) ### or?
    }
  }
  if(verbose) {
    close(pb)
  }  
  
  ### Aggregate data
  all_knee <- data.table::rbindlist(list_knee)
  all_picked_cells <- do.call(c,list_picked_cells)
  
  print(paste(
    "Number of picked cells:",
    length(all_picked_cells),
    "; note that this should be more than the actual cell count, but small enough for fast preprocessing. Final droplet selection is best done during postprocessing"))
  
  #Return result
  toret <- list(
    meta=meta,
    bascetRoot=bascetRoot,
    inputFiles=inputFiles,
    list_group_bc_count=list_group_bc_count,
    list_picked_cells=list_picked_cells,
    all_knee=all_knee,
    group_prefix=group_prefix,
    numgroup=numgroup
  )
  class(toret) <- "debstat"
  toret
}

is.debstat <- function(debstat) {
  class(debstat) == "debstat"
}



###############################################
#' Produce summary kneeplot given debarcoded statistics
#' 
#' @param debstat Statistics produced about debarcoded files
#' @param filename Optional file to store kneeplot in
#' 
#' @return A ggplot object
#' @export
DebarcodedKneePlot <- function( 
    debstat,
    filename=NULL
){
  #check arguments
  stopifnot(is.debstat(debstat))
  
  #Set line aesthetics: if kept or not
  debstat$all_knee$picked <- factor(as.character(debstat$all_knee$picked), levels = c("TRUE","FALSE"))
  
  #Set line aesthetics: input group (prefix for cell IDs)
  debstat$all_knee$prefix <- debstat$group_prefix[debstat$all_knee$group]
  
  p <- ggplot2::ggplot(debstat$all_knee, ggplot2::aes(index,count, color=prefix, linetype=picked)) + 
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10()
  
  if(!is.null(filename)) {
    ggplot2::ggsave(filename, plot=p, width = 8, height = 5)
  }
  p
}
#  DebarcodedKneePlot(debstat, filename = "kneeplot.pdf")







###############################################
#' Take debarcoded reads, merge them, and split them into suitable numbers of shards.
#' 
#' The reads from one cell is guaranteed to only be present in a single shard.
#' This makes parallel processing simple as each shard can be processed on
#' a separate computer. Using more shards means that more computers can process the data in parallel.
#' However, if you perform all the calculations on a single computer, having more
#' than one shard will not result in a speedup. This option is only relevant
#' when using a cluster of compute nodes.
#' 
#' @param debstat Plan for sharding provided by PrepareSharding
#' @param numOutputShards How many shards to generate /for each input prefix/
#' @param outputName Name of the output file: Properly sharded debarcoded reads
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @param numThreads Number of threads to use per job. Default is the number from the runner
#' @param numWriterThreads Advanced settings: Number of writer threads to use per job
#' 
#' @param totalMem How much memory to use. Extracted from runner if set
#' @param streamArenaMem Advanced settings: How much memory to use for streaming arena (fraction, given as e.g. "10%")
#' 
#' @return A runner job (details depends on runner)
#' @export
BascetShardify <- function(
    debstat,
    numOutputShards=1,  ### should it be per input group rather
    outputName="filtered", 
    overwrite=FALSE,

    numThreads=NULL,
    numWriterThreads=NULL,
    totalMem=NULL,
    streamArenaMem=NULL,
    streamBufferSize=NULL,

    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  
  #Set number of threads if not given
  if(is.null(numThreads)) {
    numThreads <- as.integer(runner@ncpu)
  }
  if(numThreads < 3) {
    print("Note: Setting number of threads to 3 as this is the minimum, even if the CPU has fewer cores")
    numThreads <- 3
  }
  
  #Check arguments 
  stopifnot(is.debstat(debstat))
  stopifnot(is.positive.integer(numOutputShards))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  streamArenaMem <- parse_size_string(streamArenaMem)
  stopifnot(is.null(streamBufferSize) || is.percent.string(streamBufferSize))
  
  
  #Figure out memory usage
  if(!is.null(totalMem)) {
    totalMem <- parse_size_string(totalMem)
    stopifnot(totalMem > fs::fs_bytes("1Gb"))
  } else {
    #Take memory from runner if possible
    if(runner@mem!="") {
      totalMem <- parse_size_string(runner@mem) - fs::fs_bytes(bascetInstance@containerMem)
      stopifnot(totalMem > fs::fs_bytes("1Gb"))
    } else {
      print("Warning: Total memory was not specified. We strongly encourage doing this to ensure performance")
    }
  }
  
  
  #Figure out mapping input vs output shards
  totalNumOutputShards <- numOutputShards*debstat$numgroup
  dfListOutputs <- data.frame(
    group=rep(1:debstat$numgroup, numOutputShards),
    outputFile=makeOutputShardNames(bascetRoot, outputName, "tirp.gz", totalNumOutputShards)
  )

  dfListInputs <- data.frame(
    group=debstat$meta$group,
    inputFile=debstat$inputFiles
  )

  dfListIO <- merge(dfListOutputs, dfListInputs)

  #Turn I/O-mapping into one list per job
  inputFiles <- list()
  outputFiles <- list()
  for(g in 1:debstat$numgroup){
    inputFiles[[g]] <- shellscriptMakeCommalist(unique(dfListIO$inputFile[dfListIO$group==g]))
    outputFiles[[g]] <- shellscriptMakeCommalist(unique(dfListIO$outputFile[dfListIO$group==g]))
  }

  all_outputFiles <- unique(dfListIO$outputFile)

  if(bascetCheckOverwriteOutput(all_outputFiles, overwrite)) {
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Zshardify",
      bascetInstance = bascetInstance,
      cmd = c(
        shellscriptMakeBashArray("files_in", inputFiles),
        shellscriptMakeBashArray("files_out", outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),
        
        shellscriptMakeFilesExpander("CELLFILE", debstat$list_picked_cells),
        assembleBascetCommand(bascetInstance, c(
          "shardify", 
          "--temp=$BASCET_TEMPDIR",
          "-i=${files_in[$TASK_ID]}",   
          "-o=${files_out[$TASK_ID]}",  
          "--include=${CELLFILE[$TASK_ID]}",

          if(!is.null(numThreads)) paste0("--threads=",numThreads),
          if(!is.null(numWriterThreads)) paste0("--numof-threads-write=",format_size_bascet(numWriterThreads)),
          if(!is.null(totalMem)) paste0("--memory=",format_size_bascet(totalMem)),
          if(!is.null(streamArenaMem)) paste0("--sizeof-stream-arena=",format_size_bascet(streamArenaMem))
        ))
      ), 
      arraysize = debstat$numgroup
    )
  } else {
    new_no_job()
  }
}







if(FALSE){
  
  #~/mystore/dataset/250611_scinfluenza/merge
  bascetRoot <- "/home/m/mahogny/mystore/dataset/250611_scinfluenza/bascet"
  
  debstat <- PrepareSharding(
    bascetRoot,
    inputName="debarcoded",
    bascetInstance=bascetInstance.default,
    minQuantile=0.99
  )  
  
  DebarcodedKneePlot(debstat, filename = "kneeplot.pdf")


  BascetShardify2(
    debstat = debstat,
    numOutputShards = 1,
    runner=SlurmRunner(bascetRunner.default, ncpu="4")  #not much CPU needed. increased for memory demands
  )
}
