
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
    minQuantile=0.9,
    bascetInstance=GetDefaultBascetInstance(),
    verbose=TRUE
){
  
  #Read metadata to know which files need to be concatenated.
  #This speeds up the process
  fstats <- file.path(bascetRoot, paste0(inputName,".meta"))
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
  inputFiles <- detectShardsForFile(bascetRoot, inputName)
  inputFiles <- naturalsort::naturalsort(inputFiles)
  
  if(length(inputFiles)!=nrow(meta)){
    stop("Metadata has different number of files than currently detected")
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
      hist_f,
      col_types = list(
        bc = readr::col_character(),
        cnt = readr::col_double()
      ),
      progress=FALSE
    )
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
    
    print(head(dat))

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
    
    #Produce a knee plot. Subsample to keep the speed high
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
  list(
    meta=meta,
    bascetRoot=bascetRoot,
    inputFiles=inputFiles,
    list_group_bc_count=list_group_bc_count,
    list_picked_cells=list_picked_cells,
    all_knee=all_knee,
    group_prefix=group_prefix,
    numgroup=numgroup
  )
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
  
  #should possibly show all points in kneeplot. currently cut off at low index, which is weird TODO
  
  #Set line aesthetics
  debstat$all_knee$picked <- factor(as.character(debstat$all_knee$picked), levels = c("TRUE","FALSE"))
  
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
#' @param numOutputShards How many shards to generate /for each input library/
#' @param outputName Name of the output file: Properly sharded debarcoded reads
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A runner job (details depends on runner)
#' @export
BascetShardify <- function(
    debstat,
    numOutputShards=1,
    outputName="filtered", 
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
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
      jobname = "Z_shardify",
      bascetInstance = bascetInstance,
      cmd = c(
        shellscriptMakeBashArray("files_in", inputFiles),
        shellscriptMakeBashArray("files_out", outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),
        
        shellscriptMakeFilesExpander("CELLFILE", debstat$list_picked_cells),
        paste(
          bascetInstance@prependCmd,
          bascetInstance@bin,
          "shardify", 
          ### "-t $BASCET_TEMPDIR", #no longer part
          ### "-@ numthreads",  #should autodetect
          "-i ${files_in[$TASK_ID]}",   
          "-o ${files_out[$TASK_ID]}",  
          "--include ${CELLFILE[$TASK_ID]}",
          "--buffer-size 16000",
          "--page-size 32"
        )  
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



