

###############################################
#' Function template, where basic parameter documentation can be obtained from
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance Configuration for how to run the Bascet Rust API
#' @param num_output_shards Number of output shards, i.e., how many output files to split the data into? (>=1)
#' @param includeCells List of cells to process
#' @param bascetFile Handle for an opened Bascet file
#' @param cellID ID of the cell (string)
#' @param verbose description
#' 
#' @return A job to be executed
template_BascetFunction <- function(
    bascetRoot, 
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()){}



###############################################
#' A wrapper to cache a computation. Put your function in as an argument,
#' as R will only compute its value if needed. If the cache file exist,
#' it will not be run again
#' 
#' @param fname Name of the file to store the cache in. The extension .RDS is added automatically
#' @return The value
BascetCacheComputation <- function(bascetRoot, fname, value){
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
#' @export
#' @return A data frame with metadata for the raw input files
DetectRawFileMeta <- function(rawRoot, verbose=FALSE){
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
#' Generate BAM with barcodes from input raw FASTQ
#' 
#' @inheritParams template_BascetFunction
#' @param rawmeta Metadata for the raw FASTQ input files. See DetectRawFileMeta
#' @param outputName Name output files: Debarcoded reads
#' @param outputNameIncomplete Name of output files: Reads that could not be parsed
#' @param chemistry The type of data to be parsed
#' @param barcodeTolerance Optional: Number of mismatches allowed in the barcode for it to still be considered valid
#' @param overwrite description
#' @export
BascetGetRaw <- function(
    bascetRoot, 
    rawmeta, 
    outputName="debarcoded", 
    outputNameIncomplete="incomplete_reads", 
    chemistry="atrandi_wgs",  #or atrandi_rnaseq; any way to get list from software?
    barcodeTolerance=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){

  #One output per input pair of reads  
  num_shards <- nrow(rawmeta)
  if(num_shards==0){
    stop("No input files")
  }
  outputFilesComplete <- makeOutputShardNames(bascetRoot, outputName, "tirp.gz", num_shards)
  outputFilesIncomplete <- makeOutputShardNames(bascetRoot, outputNameIncomplete, "tirp.gz", num_shards)
  
  #Check if libnames should be added
  add_libnames <- any(rawmeta$prefix!="")
  
  #Write a file describing the libraries
  rawmeta$shard <- 1:num_shards
  write.csv(
    rawmeta,
    file=file.path(bascetRoot, paste0(outputName, ".meta")),
    row.names = FALSE
  )

  if(bascetCheckOverwriteOutput(outputFilesComplete, overwrite)) {
    RunJob(
      runner = runner, 
      jobname = "Z_getraw",
      bascetInstance = bascetInstance,
      cmd = c(
        #shellscript_set_tempdir(bascetInstance),
        shellscriptMakeBashArray("files_r1",file.path(rawmeta$dir, rawmeta$r1)),
        shellscriptMakeBashArray("files_r2",file.path(rawmeta$dir, rawmeta$r2)),
        shellscriptMakeBashArray("libnames",rawmeta$prefix),
        shellscriptMakeBashArray("files_out",outputFilesComplete),
        shellscriptMakeBashArray("files_out_incomplete",outputFilesIncomplete),
        
        ### Abort early if needed    
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),
        
        paste(
          bascetInstance@prependCmd,
          bascetInstance@bin, 
          "getraw",
          "-t $BASCET_TEMPDIR",
          "--chemistry",chemistry,  
          if(!is.null(barcodeTolerance)) c("barcode-tol", barcodeTolerance),
          "--r1 ${files_r1[$TASK_ID]}",  
          "--r2 ${files_r2[$TASK_ID]}",
          if(add_libnames) "--libname ${libnames[$TASK_ID]}",
          "--out-complete   ${files_out[$TASK_ID]}",                 #Each job produces a single output
          "--out-incomplete ${files_out_incomplete[$TASK_ID]}"                 #Each job produces a single output
        )
      ),
      arraysize = nrow(rawmeta)
    )    
  } else {
    new_no_job()
  }
}





################################################################################
################ Bascet command line tools: WGS ################################
################################################################################




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
#' 
#' @export
BascetShardifyOld <- function(
    bascetRoot, 
    inputName="debarcoded", 
    includeCells=NULL, ############# TODO: get rid of this parameter; only support direct merging
    num_output_shards=1,
    outputName="filtered", 
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){

  input_shards <- detectShardsForFile(bascetRoot, inputName)
  if(length(input_shards)==0){
    stop("Found no input files")
  }

  #Include all cells if nothing else provided
  if(is.null(includeCells)){
    includeCells <- BascetCellNames(bascetRoot, inputName, bascetInstance=bascetInstance)
    includeCells <- unique(includeCells$cell) #when shardifying, we expect cells to appear more than once  -- could warn for other commands!
    print(paste("Including all the",length(includeCells), "cells"))
  }

  #Figure out which cell goes into which shard
  list_cell_for_shard <- shellscriptSplitArrayIntoListRandomly(includeCells, num_output_shards)
  
  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, input_shards)
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "tirp.gz", num_output_shards)

  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Z_shardify",
      bascetInstance = bascetInstance,
      cmd = c(
        #shellscript_set_tempdir(bascetInstance),
        shellscriptMakeBashArray("files_out", outputFiles),
        
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




###############################################
#' Transform data
#' 
#' This command enables
#' * subsetting to a list of cells
#' * converting between file formats
#' * merging shards
#' * dividing shards
#' 
#' @param num_divide description
#' @param num_merge description
#' @param out_format description
#' 
#' @export
BascetMapTransform <- function(
    bascetRoot, 
    inputName, 
    outputName,
    num_divide=1,
    num_merge=1,
    out_format="tirp.gz", ### not really!
    includeCells=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  
  
  is.integer.overequal1 <- function(x){
    return(x>=1 & x==as.integer(x))
  }
  
  #Verify input
  if(!(num_divide==1 || num_merge==1)){
    stop("Either divide or merge must be set to 1")
  }
  if(!is.integer.overequal1(num_divide) || !is.integer.overequal1(num_merge)){
    stop("Number of divisions and merges must be an integer >= 1")
  }
  
  #TODO check outformat

  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, detectShardsForFile(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }
  
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, out_format, num_shards)

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
        #shellscript_set_tempdir(bascetInstance),
        shellscriptMakeBashArray("files_in",inputFiles),
        shellscriptMakeBashArray("files_out",outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),
        
        if(produce_cell_list) shellscriptMakeFilesExpander("CELLFILE", list_cell_for_shard),
        paste(
          bascetInstance@prependCmd,
          bascetInstance@bin, 
          "transform",
          if(produce_cell_list) "--cells ${CELLFILE[$TASK_ID]}",
          #"-t $BASCET_TEMPDIR",  ##not supported
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






################################################################################
################ Quality control ###############################################
################################################################################



###############################################
#' 
#' Produce a kneeplot
#' 
#' @param adata description
#' @param max_species description
#' 
#' @return A ggplot object
#' @export
KneeplotPerSpecies <- function(
    adata, 
    max_species=NULL
) {
  strain_cnt <- adata@assays[[DefaultAssay(adata)]]$counts
  
  if(!is.null(max_species)){
    strain_cnt <- strain_cnt[order(rowSums(strain_cnt), decreasing = TRUE),]
    strain_cnt <- strain_cnt[1:min(nrow(strain_cnt), max_species),]
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
#' 
#' 
#' 
#' @param adata A Seurat object with the DefaultAssay set to species count
#' @return a ggarranged set of ggplots
#' @export
BarnyardPlotMatrix <- function(
    adata
){
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
#' @param useKrakenDB Path to KRAKEN2 database
#' @export
#' 
BascetRunFASTP <- function(
    bascetRoot,
    numLocalThreads=1,
    inputName="asfq", ######### should be able to take filtered and pipe to if needed  "filtered"  TODO
    outputName="fastp",
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  
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
          "--thread", numLocalThreads,
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




ReadBascetCountMatrix_one <- function(
    fname,
    verbose=FALSE
){
  
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
#' @return Count matrix as sparseMatrix
#' @export
ReadBascetCountMatrix <- function(
    bascetRoot, 
    inputName,
    verbose=FALSE
){
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
  #allmat
}




