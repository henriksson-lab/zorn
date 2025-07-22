
###############################################
#' Prepare to shard reads by collecting statistics about
#' each barcode, and filtering out cells with few reads
#' 
#' @param minQuantile Read count-based cutoff for inclusion in final shards
#' @param bascetInstance A Bascet instance
#' @return Statistics about the debarcoded reads
#' @export
PrepareSharding <- function(
    bascetRoot, 
    inputName, 
    minQuantile=0.9,
    bascetInstance=GetDefaultBascetInstance(),
    verbose=TRUE
){
  
  #Read metadata to know which files need to be concatenated.
  #This speeds up the process
  fstats <- file.path(bascetRoot, paste0(inputName,".meta"))
  meta <- read.csv(fstats)
  meta$group <- as.integer(factor(meta$prefix))
  
  #Get all the TIRPs, sum up the reads  
  inputFiles <- detectShardsForFile(bascetRoot, inputName)
  inputFiles <- naturalsort::naturalsort(inputFiles)
  
  if(length(inputFiles)!=nrow(meta)){
    stop("Metadata has different number of files than currently detected")
  }
  
  if(verbose) {
    print("Loading barcode counts")
  }
  
  #Read BC info for each file
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
  }
  
  if(verbose) {
    print("Merging barcode lists for each library")
  }
  list_group_bc_count <- list()
  for(g in unique(meta$group)) {
    if(verbose) {
      print(g)
    }
    
    #Gather counts for cells across files for each library (e.g. multiple lanes)
    dat <- data.table::rbindlist(list_hist[meta$shard[meta$group==g]])
    colnames(dat) <- c("cellid","onecount")
    
    #Add up cells across barcodes
    sumdat <- dat[ ,list(count=sum(onecount)), by=cellid]
    
    list_group_bc_count[[g]] <- sumdat[order(count,decreasing=TRUE)]
  }
  
  if(verbose) {
    print("Kneeplot analysis for each library")
  }
  list_knee <- list()
  list_picked_cells <- list()
  for(g in unique(meta$group)) {
    if(verbose) {
      print(g)
    }
    
    #Produce a knee plot. Subsample to keep the speed high
    knee <- list_group_bc_count[[g]]
    knee$index <- 1:nrow(knee)
    
    sub_knee <- knee[sample(1:nrow(knee), min(nrow(knee),10000)),, drop=FALSE]
    
    sub_knee$cs <- cumsum(sub_knee$count)
    sub_knee$picked <- sub_knee$cs > sum(sub_knee$count)*minQuantile
    sub_knee$group <- g
    
    #min_quantile <- 0.99
    read_cutoff <- quantile(sub_knee$count, min_quantile)
    sub_knee$picked <- sub_knee$count > read_cutoff
    
    list_knee[[g]] <- sub_knee
    list_picked_cells[[g]] <- knee$cellid[knee$count > read_cutoff]
  }
  
  #Aggregate data
  all_knee <- data.table::rbindlist(list_knee)
  all_knee$group <- as.character(all_knee$group)
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
    all_knee=all_knee
  )
}


if(FALSE){
  
  #~/mystore/dataset/250611_scinfluenza/merge
  bascetRoot <- "/home/m/mahogny/mystore/dataset/250611_scinfluenza/bascet"
  
  debstat <- PrepareSharding(
    bascetRoot,
    inputName="debarcoded",
    bascetInstance=bascetInstance.default,
    min_quantile=0.99
  )  
}




###############################################
#' Produce summary kneeplot given debarcoded statistics
#' 
#' @param debstat Statistics produced about debarcoded files
#' @param filename Optional file to store kneeplot in
#' @return A ggplot object
#' @export
DebarcodedKneePlot <- function( 
    debstat,
    filename=NULL
){
  
  #Set line aesthetics
  debstat$all_knee$picked <- factor(as.character(debstat$all_knee$picked), levels = c("TRUE","FALSE"))
  
  p <- ggplot2::ggplot(debstat$all_knee, ggplot2::aes(index,count, color=group, linetype=picked)) + 
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10()
  
  if(!is.null(filename)) {
    ggplot2::ggsave(filename, plot=p, width = 8, height = 5)
  }
  p
}


if(FALSE){
  DebarcodedKneePlot(bcstat, filename = "kneeplot.pdf")
}






###############################################
#' Take debarcoded reads and split them into suitable numbers of shards.
#' 
#' The reads from one cell is guaranteed to only be present in a single shard.
#' This makes parallel processing simple as each shard can be processed on
#' a separate computer. Using more shards means that more computers can be used.
#' 
#' If you perform all the calculations on a single computer, having more
#' than one shard will not result in a speedup. This option is only relevant
#' when using a cluster of compute nodes.
#' 
#' 
#' 
#' TODO if we have multiple input samples, is there a way to group them?
#' otherwise we will be reading more input files than needed. that said,
#' if we got an index, so if list of cells specified, it is possible to quickly figure out
#' out if a file is needed at all for a merge
#' 
#' TODO Figuring out if a file is needed can be done at "planning" (Zorn) stage
#' 
#' TODO seems faster to have a single merger that writes multiple output files if
#' cell list is not provided. if the overhead is accepted then read all input files and
#' discard cells on the fly
#' 
#' @param inputName Name of input file: Debarcoded reads
#' @param outputName Name of the output file: Properly sharded debarcoded reads
#' @param numOutputShards How many shards to generate /for each input library/
#' 
#' @export
BascetShardify2 <- function(
    useStats,
    numOutputShards=1, ### for each input library
    outputName="filtered", 
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  
  
  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  if(length(input_shards)==0){
    stop("Found no input files")
  }
  
  #Figure out which cell goes into which shard
  list_cell_for_shard <- shellscriptSplitArrayIntoListRandomly(includeCells, numOutputShards)
  
  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, input_shards)
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "tirp.gz", numOutputShards)
  
  
  
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Z_shardify",
      bascetInstance = bascetInstance,
      cmd = c(
        #shellscript_set_tempdir(bascetInstance),
        shellscript_make_bash_array("files_out", outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),
        
        shellscriptMakeFilesExpander("CELLFILE", list_cell_for_shard),
        paste(
          bascetInstance@prependCmd,
          bascetInstance@bin,
          "shardify", 
          "-t $BASCET_TEMPDIR",
          "-i",shellscriptMakeCommalist(inputFiles), #Need to give all input files for each job
          "-o ${files_out[$TASK_ID]}",                 #Each job produces a single output
          "--cells $CELLFILE"                          #Each job takes its own list of cells
        ),  
        "rm $CELLFILE"
      ), 
      arraysize = num_output_shards
    )
  } else {
    new_no_job()
  }
}

























if(FALSE) {
  BascetShardify2(
    bascetRoot,
    useStats = bcstat,
    num_output_shards = 20,
    runner=SlurmRunner(bascet_runner.default, ncpu="4")  #not much CPU needed. increased for memory demands
  )
}



