


################################################################################
################ Helper functions: memory size #################################
################################################################################



if(FALSE) {
  # new!
  # to be compatible with Rust bytesize: https://docs.rs/bytesize/latest/bytesize/
  fs::fs_bytes("3")
  fs::fs_bytes("3mb")
  ?fs::fs_bytes
  fs::fs_bytes("3g") ### need to support 
  
  as.character(fs::fs_bytes("3GB"))
  as.character(fs::fs_bytes("4Kib"))
  as.character(fs::fs_bytes("4"))
  as.character(fs::fs_bytes("4b"))
  
  inp <- "5Mb"
  
  inp_conv <- fs::fs_bytes(inp)
  if(is.na(inp_conv)) {
    #If R library fails to convert bytes, we will attempt some additional
    #patters used by SLURM
  }
  
  sorig <- "5g"
  regmatches(sorig,regexpr("[0-9]+[gGmMkKtT]$",sorig))
  regexpr("[0-9]+[gGmMkKtT]",sorig)
  
  
}



###############################################
#' Parse a string with a size, such as 1g, 1m, 1k, or just 123 (bytes)
#' 
#' @param s Size as string. If numeric, it is just returned
#' 
#' @return Size in bytes, as integer in a string (long format)
parse_size_string <- function(s) {
  if(is.null(s)) {
    return(NULL)
  }
  
  s_fs <- fs::fs_bytes(s)
  if(!is.na(s_fs)) {
    #Return fs_bytes response if successful
    return(s_fs)
  } else {
    #If R library fails to convert bytes, we will attempt some additional patterns used by SLURM.
    #These include lower case units such as "5g"
    sorig <- s
    
    #If there is a trailing b, then remove it
    if(stringr::str_ends(s, stringr::fixed("b"))) {
      s <- stringr::str_sub(s, 1,stringr::str_length(s)-1)
    }
    
    #Figure out prefix (g, etc)
    pref <- stringr::str_sub(s, 1,stringr::str_length(s)-1)
    last_c <- stringr::str_sub(s, stringr::str_length(s))

    if(last_c=="t") {
      return(fs::fs_bytes(paste0(pref,"T")))
    } else if(last_c=="g") {
      return(fs::fs_bytes(paste0(pref,"G")))
    } else if(last_c=="m") {
      return(fs::fs_bytes(paste0(pref,"M")))
    } else if(last_c=="k") {
      return(fs::fs_bytes(paste0(pref,"k")))
    } else {
      stop(paste("Cannot parse memory size:", s))
      #return(NA)
    }
#    > class(inp_conv)
#    [1] "fs_bytes" "numeric" 
  }
}
#parse_size_to_bytes("2dasd")


###############################################
#' Format size for input to Bascet
#' 
#' @param s Size as output from fs::fs_bytes, or a string
format_size_bascet <- function(s) {
  paste0(s, "B")
}

###############################################
#' Parse a string with a size, such as 1g, 1m, 1k, or just 123 (bytes)
#' 
#' @param s Size as string. If numeric, it is just returned
#' 
#' @return Size in mb, as integer in a string (long format)
#parse_size_to_rust_mb <- function(s) {
  #print(777)
  #print(s)
#  b <- round(parse_size_to_bytes(s)/1000000)
#  b
  #  format(b, scientific =FALSE)
#}

formatPlainNumber <- function(s) {
  format(s, scientific=FALSE)
}


###############################################
#' Given memory amount in mb, format it in a format suitable for Rust parsing
#' 
#formatMemMB <- function(s) {
#  paste0(format(s, scientific=FALSE),"mb")
#}


###############################################
#' Helper for functions. Set total memory argument based on container and
#' runner, if not user specified
#' 
#' @param totalMem Total memory to allocate, as a string (e.g. "8g"). If NULL, derived from runner
#' @param runner The job runner, used to extract memory settings if totalMem is NULL
#' @param bascetInstance A Bascet instance, used to subtract container memory overhead
checkTotalMemArg <- function(
    totalMem,
    runner,
    bascetInstance
) {
  #Check memory sizes
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
  totalMem
}


################################################################################
################ Helper functions: check if call is correct ####################
################################################################################

###############################################
#' Check that parameter is in the form "xxx%"
#' @param s A string to test for percent format (e.g. "10%")
is.percent.string <- function(s) {
  if(!is.null(s)) {
    pref <- stringr::str_sub(s, 1,stringr::str_length(s)-1)
    last_c <- stringr::str_sub(s, stringr::str_length(s))
    
    if(last_c=="%" & is.numeric(pref)) {
      val <- as.numeric(pref)
      val >= 0 && val<=100
    } else {
      FALSE
    }
  }
}


###############################################
#' Check that parameter is a valid memory size
#' @param x A string representing a memory size (e.g. "8g")
is.valid.memsize <- function(x) {
  !is.na(parse_size_string(x))
}


###############################################
#' Check that parameter is a valid thread count
#' @param x A numeric thread count
is.valid.threadcount <- function(x) {
  #Note: not calling is.positive.integer to ensure we get a proper error message
  round(x)==x & x>0
}

###############################################
#' Check that parameter is an integer and >0
#' @param x A numeric value
is.positive.integer <- function(x) {
  round(x)==x & x>0
}

###############################################
#' Check that parameter is castable to an integer
#' @param x A numeric value
is.integer.like <- function(x) {
  round(x)==x
}

###############################################
#' Check that parameter is a valid shard name
#' @param x A string representing a shard name
is.valid.shardname <- function(x) {
  # can expand upon this
  is.character(x) && !stringr::str_detect(x,stringr::fixed("."))
}


###############################################
#' Check that parameter is a valid list of cells
#' @param x A character vector of cell names, or NULL
is.valid.listcells <- function(x) {
  is.null(x) || is.character(x)
}


###############################################
#' Check that parameter is a valid shard name
#' @param x A file path string pointing to a FASTA file
is.existing.fasta <- function(x) {
  is_fasta <- 
    stringr::str_ends(x,stringr::fixed(".fa")) ||
    stringr::str_ends(x,stringr::fixed(".fasta")) ||
    stringr::str_ends(x,stringr::fixed(".fa.gz")) ||
    stringr::str_ends(x,stringr::fixed(".fasta.gz"))
  file.exists(x) && is_fasta
}


###############################################
#' Check that parameter is a number between 0..1
#' @param x A numeric value
is.numeric.range01 <- function(x) {
  is.numeric(x) && x>=0 && x<=1
}





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
#' @param outputFiles Files that we expect to exist
#' @param overwrite Files that we expect to exist
#' 
#' @return TRUE if ok to proceed
bascetCheckOverwriteOutput <- function(
  outputFiles, 
  overwrite
){
  if(all(file.exists(outputFiles)) & !overwrite){
    print("All files to write already exist; skipping. To change this behaviour, set overwrite=TRUE")
    FALSE
  } else {
    TRUE
  }
}

  

################################################################################
################ Bascet caching system #########################################
################################################################################




###############################################
#' A wrapper to cache a computation. Put your function in as an argument,
#' as R will only compute its value if needed. If the cache file exist,
#' it will not be run again
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param fname Name of the file to store the cache in. The extension .RDS is added automatically
#' @param value The value to be cached
#' 
#' @return The cached value
#' @export
BascetCacheComputation <- function(
    bascetRoot, 
    fname, 
    value
){
  #Check input arguments 
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.character(fname))

  
  fname <- file.path(bascetRoot,paste0(fname,".RDS"))
  if(file.exists(fname)){
    print("Found previously cached value")
    readRDS(fname)
  } else {
    print("Running calculation and caching value")
    saveRDS(value, fname) 
    value
  }
}



###############################################
#' Transform data
#' 
#' This command enables
#' * subsetting to a list of cells
#' * converting between file formats
#' * merging shards
#' * dividing shards
#' 
#' 
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param numDivide Must be >=1; if more than 1, divide each input shard into number of outputs given here
#' @param numMerge Must be >=1; if more than 1, merge this number of input shards into one
#' @param outFormat Extension for the output files
#' @param includeCells List of cells to include, or NULL if to include all
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetMapTransform <- function(
    bascetRoot, 
    inputName, 
    outputName,
    numDivide=1,
    numMerge=1,
    outFormat="tirp.gz", ### not really!
    includeCells=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  #Check input arguments 
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.positive.integer(numDivide))
  stopifnot(is.positive.integer(numMerge))
  if(!(numDivide==1 || numMerge==1)){
    stop("Either divide or merge must be set to 1")
  }
  stopifnot(is.valid.listcells(includeCells))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  
  #TODO check outformat

  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, detectShardsForFile(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }
  
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, outFormat, num_shards)

  #If cell list is provided, produce a file for input (not all transform calls can handle this, so optional)
  produce_cell_list <- !is.null(includeCells)
  if(produce_cell_list) {
    #Figure out which cell goes into which shard   TODO, using the same cell list for all shards
    list_cell_for_shard <- list()
    for(i in 1:length(inputFiles)){
      list_cell_for_shard[[i]] <- includeCells
    }
  }
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    #Run the job
    RunJob(
      runner = runner, 
      jobname = "Ztransform",
      bascetInstance = bascetInstance,
      cmd = c(
        shellscriptMakeBashArray("files_in",inputFiles),
        shellscriptMakeBashArray("files_out",outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),
        
        if(produce_cell_list) shellscriptMakeFilesExpander("CELLFILE", list_cell_for_shard),
        assembleBascetCommand(bascetInstance, c(
          "transform",
          if(produce_cell_list) "--cells=${CELLFILE[$TASK_ID]}",
          "-i=${files_in[$TASK_ID]}",
          "-o=${files_out[$TASK_ID]}"
        ))
      ),
      arraysize = num_shards
    )  
  } else {
    new_no_job()
  }
}




###############################################
#' Convert data to Bascet-FASTQ
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param numLocalThreads Number of threads to use per job. Default is the number from the runner
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetToFastq <- function(
    bascetRoot, 
    inputName, 
    outputName,
    numLocalThreads=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  
  #Set number of threads if not given
  if(is.null(numLocalThreads)) {
    numLocalThreads <- as.integer(runner@ncpu)
  }
  stopifnot(is.valid.threadcount(numLocalThreads))
  
  #Check input arguments 
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  
  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, detectShardsForFile(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }
  
  outputFilesR1 <- makeOutputShardNames(bascetRoot, outputName, "R1.fq.gz", num_shards)
  outputFilesR2 <- makeOutputShardNames(bascetRoot, outputName, "R2.fq.gz", num_shards)

  if(bascetCheckOverwriteOutput(c(outputFilesR1, outputFilesR2), overwrite)) {
    #Run the job
    RunJob(
      runner = runner, 
      jobname = "Ztofq",
      bascetInstance = bascetInstance,
      cmd = c(
        shellscriptMakeBashArray("files_in",inputFiles),
        shellscriptMakeBashArray("files_out_r1",outputFilesR1),
        shellscriptMakeBashArray("files_out_r2",outputFilesR2),
        
        ### Abort early if needed    
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out_r1[$TASK_ID]}"),
        
        assembleBascetCommand(bascetInstance, c(
          "to-fastq",
          "-i=${files_in[$TASK_ID]}",
          "--out-r1=${files_out_r1[$TASK_ID]}",
          "--out-r2=${files_out_r2[$TASK_ID]}",
          "--temp=$BASCET_TEMPDIR",
          paste0("-@=",numLocalThreads) 
        ))
      ),
      arraysize = num_shards
    )  
  } else {
    new_no_job()
  }
}





################################################################################
################ Quality control ###############################################
################################################################################



###############################################
#' 
#' Produce a kneeplot
#' 
#' @param adata A Seurat object with the DefaultAssay having counts per species (or similar)
#' @param maxSpecies Maximum number of species to show. The most abundant species will be shown first
#' 
#' @return A ggplot object
#' @export
KneeplotPerSpecies <- function(
    adata, 
    maxSpecies=NULL
) {
  #check arguments
  #TODO

  strain_cnt <- adata@assays[[DefaultAssay(adata)]]$counts
  
  if(!is.null(maxSpecies)){
    strain_cnt <- strain_cnt[order(rowSums(strain_cnt), decreasing = TRUE),]
    strain_cnt <- strain_cnt[1:min(nrow(strain_cnt), maxSpecies),]
  }
  
  allknee <- list()
  for(i in 1:nrow(strain_cnt)){
    onedf <- data.frame(
      strain = rownames(strain_cnt)[i],
      cnt = strain_cnt[i,]
    )
    onedf <- onedf[order(onedf$cnt, decreasing = TRUE),]
    onedf$index <- 1:nrow(onedf)
    allknee[[paste("s",i)]] <- onedf
  }
  allknee <- do.call(rbind, allknee)
  
  ggplot(allknee, aes(index, cnt, color=strain)) + geom_line() +
    scale_x_log10() +
    scale_y_log10() +
    xlab("Log10 Cell index") +
    ylab("Log10 Total read count") +
    theme_bw()
  # +
  #  theme(legend.position = "none")
  
}




###############################################
#' 
#' Produce a matrix of Barnyard plots, i.e., counts for one species vs another, 
#' for all combinations of species.
#' 
#' @param adata A Seurat object with the DefaultAssay holding counts per species (or similar)
#' 
#' @return a ggarranged set of ggplots
#' @export
BarnyardPlotMatrix <- function(
    adata
){
  #check arguments
  #TODO
  
  
  cnt <- adata@assays[[DefaultAssay(adata)]]$counts
  cnt <- cnt[rowSums(cnt)>0,]  ##Only consider species we have
  list_species <- rownames(cnt)#[1:5] ######## Just do a few!
  
  all_plots <- list()
  for(i in seq_along(list_species)){
    for(j in seq_along(list_species)){
      #print(paste(i,j))
      if(i>=j) {
        all_plots[[paste(i,j)]] <- ggplot()
      } else {
        df <- data.frame(
          x=cnt[i,],
          y=cnt[j,]
        )
        p <- ggplot(df, aes(x+1,y+1)) +
          geom_point() + 
          scale_x_log10() + 
          scale_y_log10() + 
          theme_bw()+
          xlab(paste("Pseudocount", list_species[i]))+
          ylab(paste("Pseudocount", list_species[j]))
        
        all_plots[[paste(i,j)]] <- p        
      }
    }
  }
  egg::ggarrange(plots = all_plots, nrow = length(list_species))  
}








###############################################
#' Run FASTP for each cell. Input must be in FASTQ file format
#'
#'
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param numLocalThreads Number of threads to use for FASTP. Default is the maximum, taken from runner settings
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#'
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetRunFASTP <- function(
    bascetRoot,
    inputName="asfq", ######### should be able to take filtered and pipe to if needed  "filtered"  TODO
    outputName="fastp",
    numLocalThreads=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  #Set number of threads if not given
  if(is.null(numLocalThreads)) {
    numLocalThreads <- as.integer(runner@ncpu)
  }
  
  #check arguments
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.valid.threadcount(numLocalThreads))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  
  #Figure out input and output file names  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles_R1 <- file.path(bascetRoot, input_shards)
  
  outputFiles_R1 <- makeOutputShardNames(bascetRoot, outputName, "R1.fq.gz", num_shards) 
  outputFiles_R2 <- makeOutputShardNames(bascetRoot, outputName, "R2.fq.gz", num_shards) 

  outputFiles_report_json <- makeOutputShardNames(bascetRoot, outputName, "html", num_shards) 
  outputFiles_report_html <- makeOutputShardNames(bascetRoot, outputName, "json", num_shards) 
  
  ### Check if paired or not
  is_paired <- isPairedFastq(inputFiles_R1[1])
  print(paste("Detect paired FASTQ:",is_paired))
  
  ### Figure out R2 names
  if(is_paired){
    inputFiles_R2 <- getFastqR2fromR1(inputFiles_R1)
  }
  
  
  if(bascetCheckOverwriteOutput(outputFiles_R1, overwrite)) {
    #Run the job
    RunJob(
      runner = runner, 
      jobname = paste0("Zfastp"),
      bascetInstance = bascetInstance,
      cmd = c(
        #shellscript_set_tempdir(bascetInstance),
        shellscriptMakeBashArray("files_html",outputFiles_report_json),
        shellscriptMakeBashArray("files_json",outputFiles_report_html),
        shellscriptMakeBashArray("files_in_R1",inputFiles_R1),
        if(is_paired) shellscriptMakeBashArray("files_in_R2",inputFiles_R2),
        shellscriptMakeBashArray("files_out_R1",outputFiles_R1),
        if(is_paired) shellscriptMakeBashArray("files_out_R2",outputFiles_R2),
        
        ### Abort early if needed    
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out_R1[$TASK_ID]}"),
        
        paste(
          bascetInstance@prependCmd,
          "fastp",
          "--thread", numLocalThreads,
          "-g",                           #polyG trimming
          "-h ${files_html[$TASK_ID]}",
          "-j ${files_json[$TASK_ID]}",
          "-i ${files_in_R1[$TASK_ID]}",
          if(is_paired) "-I ${files_in_R2[$TASK_ID]}",
          "-o ${files_out_R1[$TASK_ID]}",
          if(is_paired) "-O ${files_out_R2[$TASK_ID]}"
        )
      ),
      arraysize = num_shards
    )  
  } else {
    new_no_job()
  }
}










