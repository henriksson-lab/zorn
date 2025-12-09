


################################################################################
################ Helper functions: check if call is correct ####################
################################################################################



###############################################
#' Check that parameter is a valid thread count
is.valid.threadcount <- function(x) {
  #Note: not calling is.positive.integer to ensure we get a proper error message
  round(x)==x & x>0
}

###############################################
#' Check that parameter is an integer and >0
is.positive.integer <- function(x) {
  round(x)==x & x>0
}

###############################################
#' Check that parameter is castable to an integer
is.integer.like <- function(x) {
  round(x)==x
}

###############################################
#' Check that parameter is a valid shard name
is.valid.shardname <- function(x) {
  # can expand upon this
  is.character(x) && !stringr::str_detect(x,stringr::fixed("."))
}


###############################################
#' Check that parameter is a valid list of cells
is.valid.listcells <- function(x) {
  is.null(x) || is.character(x)
}


###############################################
#' Check that parameter is a valid shard name
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
is.numeric.range01 <- function(x) {
  is.numeric(x) && x>=0 && x<=1
}



###############################################
#' Parse a string with a size, such as 1g, 1m, 1k, or just 123 (bytes)
#' 
#' @param s Size as string. If numeric, it is just returned
#' 
#' @return Size in bytes, as integer
parse_size_to_bytes <- function(s) {
  #Return a number directly
  if(is.numeric(s)) {
    return(s)
  } 
  
  #Parse string
  sorig <- s
  mult <- 1
  pref <- stringr::str_sub(s, 1,stringr::str_length(s)-1)
  last_c <- stringr::str_sub(s, stringr::str_length(s))
  #print(last_c)
  
  if(last_c=="g") {
    mult <- 1e9
    s <- pref
  } else if(last_c=="m") {
    mult <- 1e6
    s <- pref
  } else if(last_c=="b") {
    mult <- 1e3
    s <- pref
  }
  asint <- as.integer(s)
  if(is.na(asint)) {
    stop(paste("Failed to parse size:",sorig))
  }
  asint*mult
}
#parse_size_to_bytes("2dasd")


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
      jobname = "Z_transform",
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
#' @param numLocalThreads Number of threads to use for FASTP. Default is the maximum, taken from runner settings
#' @param inputName Name of input shard
#' @param outputName Name of output shard
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
      jobname = paste0("Z_fastp"),
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
          "--thread", numLocalThreads,  ## TODO should there be = here?
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







################################################################################
############################# Count matrix / anndata ###########################
################################################################################



### Internal helper function
ReadBascetCountMatrix_one <- function(
    fname,
    verbose=FALSE
){
  #check arguments
  stopifnot(file.exists(fname))
  stopifnot(is.logical(verbose))
  
  #fname <- "/husky/henriksson/atrandi//v4_wgs_novaseq1/chromcount.1.h5"
  #fname <- "/home/mahogny/test/cnt_feature.hdf5"
  #fname <- "/husky/henriksson/atrandi//v4_wgs_novaseq1/kmer_counts.1.h5"
  
  h5f <- rhdf5::H5Fopen(fname)
  indices <- h5f$X$indices + 1 
  indptr <-  h5f$X$indptr
  dat <- h5f$X$data
  shape <- h5f$X$shape
  
  #shape
  #print(indices)
  
  #print(paste0("Assembling matrix, size: ", shape[1],"x",shape[2]))
  mat <- Matrix::sparseMatrix(  
    j=indices,   #i??   was j when fine
    p=indptr,
    x=dat,
    dims=h5f$X$shape
  )
  
  rownames(mat) <- h5f$obs$`_index`  #names of cells
  colnames(mat) <- h5f$var$`_index`  #feature names
  
  ### Read obs matrix
  
  #print(666)
  #print(names(h5f$obs))
  #list_obs_col <- setdiff(names(h5f$obs),"_index")
  #print(h5f$obs["_unmapped"])
  
  #df_obs <- data.frame(row.names=h5f$obs$`_index`)
  df_obs <- as.data.frame(h5f$obs)
  colnames(df_obs) <- names(h5f$obs)
  #print(head(df_obs))
  
  rhdf5::H5close()
  
  list(
    obs=df_obs,
    X=mat
  )
  #mat
}


###############################################
#' Read a count matrix as produced by Bascet (hdf5 format).
#' This can be output from both BascetQueryFq and BascetCountChrom
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param verbose Print additional information, primarily to help troubleshooting
#' 
#' @return Count matrix as sparseMatrix
#' @export
ReadBascetCountMatrix <- function(
    bascetRoot, 
    inputName,
    verbose=FALSE
){
  #check arguments
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.logical(verbose))

  print("Loading HDF5 files")
  
  #Figure out input file names  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards)
  if(tools::file_ext(inputFiles[1])!="h5"){
    stop("Wrong input format. should be hd5")
  }

  #Show a progress bar  
  pbar <- progress::progress_bar$new(total = length(inputFiles)*2)
  pbar$tick(0)
  
  #Load individual matrices. Sizes may not match
  list_mat <- list()
  list_obs <- list()
  for(f in inputFiles){
    list_one <- ReadBascetCountMatrix_one(f)
    
    if(verbose){
      print(dim(list_one$X))
    }
    list_mat[[f]] <- list_one$X
    list_obs[[f]] <- list_one$obs
    pbar$tick()
  }
  
  #Find union of features  
  all_colnames <- sort(unique(unlist(lapply(list_mat, colnames))))
  if(verbose){
    print(all_colnames)
  }
  num_col <- length(all_colnames)
  map_name_to_i <- data.frame(row.names = all_colnames, ind=1:length(all_colnames))
  if(verbose){
    print(map_name_to_i)
  }
  
  #Make sizes compatible
  list_resized_mat <- list()
  for(f in inputFiles){
    mat <- list_mat[[f]]
    new_mat <- MatrixExtra::emptySparse(nrow = nrow(mat), ncol = num_col, format = "R", dtype = "d")
    new_mat[1:nrow(mat), map_name_to_i[colnames(mat),]] <- MatrixExtra::as.csr.matrix(mat)  #manually look up column names!  #here, x[.,.] <- val : x being coerced from Tsparse* to CsparseMatrix
    rownames(new_mat) <- rownames(mat)
    colnames(new_mat) <- all_colnames
    # print(dim(new_mat))
    list_resized_mat[[f]] <- new_mat
    pbar$tick()
  }
  
  #Concatenate matrices
  allmat <- do.call(rbind, list_resized_mat) #TODO check that above worked properly!
  
  #Concat obs
  allobs <- do.call(rbind, list_obs)
 
  list(
    X=allmat,
    obs=allobs
  ) 
}

