



###############################################
#' Select kmers that appear useful for clustering
#' 
#' 
#' CHANGED: it just sums up kmers
#' 
#' @inheritParams template_BascetFunction
#' @return TODO
#' @export
BascetFeaturiseKMC <- function( ########### need a better name; KMC something?
    bascetRoot, 
    inputName="kmc", 
    outputName="sumkmers", 
    includeCells=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascet_instance=GetDefaultBascetInstance()
){
  
  
  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, detect_shards_for_file(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }

  outputFiles <- make_output_shard_names(bascetRoot, outputName, "kmcdb", num_shards)
  
  #If cell list is provided, produce a file for input (not all transform calls can handle this, so optional)
  produce_cell_list <- !is.null(includeCells)
  if(produce_cell_list) {
    #Currently using the same cell list for all shards  (good idea?)
    list_cell_for_shard <- list()
    for(i in 1:length(inputFiles)){
      list_cell_for_shard[[i]] <- includeCells
    }
  }
  
  if(bascet_check_overwrite_output(outputFiles, overwrite)) {
    #Run the job
    RunJob(
      runner = runner, 
      jobname = "bascet_query",
      cmd = c(
        shellscript_set_tempdir(bascet_instance),
        if(produce_cell_list) shellscript_make_files_expander("CELLFILE", list_cell_for_shard),
        shellscript_make_bash_array("files_in",inputFiles),
        shellscript_make_bash_array("files_out",outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) helper_cancel_job_if_file_exists("${files_out[$TASK_ID]}"),

        paste(
          bascet_instance@prepend_cmd,
          bascet_instance@bin, 
          "featurise",
          if(produce_cell_list) "--cells $CELLFILE",
          "-t $BASCET_TEMPDIR",
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
#' Build count table from kmer table and a list of selected kmers
#' 
#' @inheritParams template_BascetFunction
#' @param useKMERs description
#' @return TODO
#' @export
BascetQueryKMC <- function(
    bascetRoot, 
    inputName="kmc",
    outputName="kmer", 
    useKMERs,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascet_instance=GetDefaultBascetInstance()
){
  
  
  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, detect_shards_for_file(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }

  if(length(useKMERs)==0){
    stop("No KMERs")
  }

  outputFiles <- make_output_shard_names(bascetRoot, outputName, "counts.hdf5", num_shards)
  
  if(bascet_check_overwrite_output(outputFiles, overwrite)) {
    #Run the job
    RunJob(
      runner = runner, 
      jobname = "bascet_query",
      cmd = c(
        shellscript_set_tempdir(bascet_instance),
        shellscript_make_bash_array("files_in",inputFiles),
        shellscript_make_bash_array("files_out",outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) helper_cancel_job_if_file_exists("${files_out[$TASK_ID]}"),

        shellscript_make_one_file_expander("KMERFILE", useKMERs), 
        paste(
          bascet_instance@prepend_cmd,
          bascet_instance@bin, 
          "query",
          "-t $BASCET_TEMPDIR",
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







################################################################################
########################## Count matrix tools ##################################
################################################################################


###############################################
#' Read a count matrix as produced by Bascet (hdf5 format).
#' This can be output from both BascetQueryFq and BascetCountChrom
#' 
#' @return Count matrix as sparseMatrix
#' @export
ReadBascetCountMatrix <- function(
    bascetRoot, 
    inputName
){
  print("Loading HDF5 file")
 
  
  #Figure out input file names  
  input_shards <- detect_shards_for_file(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards)
  if(tools::file_ext(inputFiles[1])!="hd5"){
    stop("Wrong input format. should be hd5")
  }
  
  #Load individual matrices. Sizes may not match
  list_mat <- list()
  for(f in inputFiles){
    mat <- ReadBascetCountMatrix_one(f)
    #print(dim(mat))
    #print(colnames(mat))
#    print(table(colnames(mat)))
    list_mat[[f]] <- mat
  }

  #Find union of features  
  all_colnames <- sort(unique(unlist(lapply(list_mat, colnames))))
  print(all_colnames)
  num_col <- length(all_colnames)
  map_name_to_i <- data.frame(row.names = all_colnames, ind=1:length(all_colnames))
  print(map_name_to_i)
  
  #Make sizes compatible
  list_resized_mat <- list()
  for(f in inputFiles){
    mat <- list_mat[[f]]
    new_mat <- MatrixExtra::emptySparse(nrow = nrow(mat), ncol = num_col, format = "R", dtype = "d")
    rownames(new_mat) <- rownames(mat)
    colnames(new_mat) <- all_colnames
    #new_mat[1:nrow(mat), colnames(mat)] <- MatrixExtra::as.csr.matrix(mat)
    #print(map_name_to_i[colnames(mat),])
    new_mat[1:nrow(mat), map_name_to_i[colnames(mat),]] <- MatrixExtra::as.csr.matrix(mat)  #manually look up column names!
   # print(dim(new_mat))
    list_resized_mat[[f]] <- new_mat
  }
  
  #Concatenate matrices
  do.call(rbind, list_resized_mat) #TODO check that above worked properly!
}



ReadBascetCountMatrix_one <- function(
    fname
){
  #print("Loading HDF5 file")
  h5f <- rhdf5::H5Fopen(fname)
  indices <- h5f$X$indices+1
  indptr <-  h5f$X$indptr
  dat <- h5f$X$data
  shape <- h5f$X$shape
  
  #print(paste0("Assembling matrix, size: ", shape[1],"x",shape[2]))
  mat <- Matrix::sparseMatrix(  
    j=indices, 
    p=indptr,
    x=dat,
    dims=h5f$X$shape
  )
  
  colnames(mat) <- h5f$var$`_index`
  rownames(mat) <- h5f$obs$`_index`
  
  rhdf5::H5close()
  
  t(mat)
}




################################################################################
############################# KMC KMER tools ###################################
################################################################################


###############################################
#' Generate a histogram from a KMC3 database
#' 
#' @param fname Full name of KMC3 database
#' @param toplot Decide if to plot or return raw data
#' @return A ggplot object if toplot=TRUE, otherwise a data.frame
#' @export
KmcGetHistogram <- function(
    fname, 
    plot=TRUE
){
  library(ggplot2)
  tfile <- tempfile()
  system(paste("kmc_tools transform ",fname," histogram", tfile))
  dat <- read.table(tfile)
  colnames(dat) <- c("num_kmer","freq")
  if(plot){
    ggplot2::ggplot(dat, ggplot2::aes(log10(1+num_kmer), log10(1+freq))) + 
      ggplot2::geom_line() + 
      ggplot2::theme_bw()
  } else {
    dat
  }
}



###############################################
#' Get KMERs from a KMC3 database
#' 
#' TODO finish properly
#' 
#' @return TODO
#' @export
KmcGetKmers <- function(
    fname, 
    mincount=NULL, 
    maxcount=NULL
) {
  
  #-ci<value> - minimum value of counter to be stored in the otput file
  #-cx<value> - maximum value of counter to be stored in the otput file
  
  tfile <- tempfile()
  system(paste(
    "kmc_tools",
    "transform",
    "/home/mahogny/github/bascet/testdata/all_kmc",
    if(!is.null(mincount)) paste0("-ci",mincount),
    if(!is.null(maxcount)) paste0("-cx",maxcount),
    "dump", 
    tfile
  ))
  dat <- read.table(tfile)
  colnames(dat) <- c("kmer","freq")
  dat
}




###############################################
#' Pick random KMERs from KMC3 database. The choice is among KMERs within a frequency range
#' 
#' @param fname description
#' @param num_pick description
#' @param minfreq description
#' @param maxfreq description
#' @return List of KMERs
#' @export
KmcChooseKmerFeatures <- function(
    fname, 
    num_pick=1000, 
    minfreq=0.01, 
    maxfreq=0.10
) {
  
  ## Possibly expensive to get them all... is the total count stored somewhere? can we sample?
  dat <- KmcGetKmers(fname)
  dat <- dat[order(dat$freq),]
  #total_count <- sum(dat$freq)/nrow(dat) ##not sure if this is the best way
  num_kmer <- nrow(dat)
  
  all_kmer <- dat$kmer[round(num_kmer*minfreq):round(num_kmer*maxfreq)]
  
  
  #  
  #mincount <- round(minfreq*total_count)
  #maxcount <- round(maxfreq*total_count)
  #all_kmer <- dat$kmer[dat$freq>mincount & dat$freq<maxcount]

  if(length(all_kmer)<num_pick){
    stop(paste("Cannot pick",num_pick," KMERS; only",length(all_kmer),"pass frequency criteria"))  
  }
    
  set.seed(0)
  sample(all_kmer, num_pick)
  
  #use_kmers <- sample(KmcGetKmers("/home/mahogny/github/bascet/testdata/all_kmc", mincount = 5, maxcount = 10)$kmer, 1000)
}
  





###############################################
###############################################
###############################################
if(FALSE){
  
  ### TODO command to run KMC for query
  ### TODO command to make count matrix
  ### TODO command to make histogram and select kmers

  ## used ??  
  CreateKmerAssay <- function(counts) {
    chrom_assay <- CreateAssayObject(  ##################### in the future, can add other metadata in here too for visualization?
      counts = counts
    )
    chrom_assay
  }

  KmcGetHistogram("/home/mahogny/github/bascet/testdata/all_kmc")
  
  KmcGetKmers("/home/mahogny/github/bascet/testdata/all_kmc")
  KmcGetKmers("/home/mahogny/github/bascet/testdata/all_kmc", mincount = 5, maxcount = 10)

  set.seed(0)
  use_kmers <- sample(KmcGetKmers("/home/mahogny/github/bascet/testdata/all_kmc", mincount = 5, maxcount = 10)$kmer, 1000)
  
  KmcChooseKmerFeatures("/home/mahogny/github/bascet/testdata/all_kmc")
}



###############################################
############################################### de novo kmer system, new
###############################################



###############################################
#' Compute minhashes for each cell.
#' This is a thin wrapper around BascetMapCell
#' 
#' @inheritParams template_BascetFunction
#' @return TODO
#' @export
BascetComputeMinhash <- function( 
    bascetRoot, 
    inputName="filtered", 
    outputName="minhash", 
    #includeCells=NULL,
    overwrite=FALSE,
    max_reads=100000,  #This is most likely enough to get an overall histogram
    kmer_size=31,
    runner=GetDefaultBascetRunner(),
    bascet_instance=GetDefaultBascetInstance()
){
  BascetMapCell(
    bascetRoot=bascetRoot, 
    withfunction="_minhash_fq", 
    inputName=inputName, 
    outputName=outputName,
    #includeCells=includeCells
    args = list(
      KMER_SIZE=format(kmer_size, scientific=FALSE),
      MAX_READS=format(max_reads, scientific=FALSE)
    ),
    overwrite=overwrite,
    runner=runner,
    bascet_instance=bascet_instance)
}




###############################################
#' Gather all minhashes into a single histogram file
#' 
#' @inheritParams template_BascetFunction
#' @return A job
#' @export
BascetMakeMinhashHistogram <- function( 
  bascetRoot, 
  inputName="minhash", 
  outputName="minhash_hist.csv", 
  includeCells=NULL,
  overwrite=FALSE,
  runner=GetDefaultBascetRunner(),
  bascet_instance=GetDefaultBascetInstance()
){

  #Figure out input and output file names
  inputFiles <- file.path(bascetRoot, detect_shards_for_file(bascetRoot, inputName))
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

  if(bascet_check_overwrite_output(outputFile, overwrite)) {
    #Make the command
    cmd <- c(
      shellscript_set_tempdir(bascet_instance),
      if(produce_cell_list) shellscript_make_files_expander("CELLFILE", list_cell_for_shard),
      paste(
        bascet_instance@prepend_cmd,
        bascet_instance@bin, 
        "minhash-hist",
        if(produce_cell_list) "--cells $CELLFILE",
        "-t $BASCET_TEMPDIR",
        "-i", shellscript_make_commalist(inputFiles),
        "-o", outputFile
      )
    )
  
    #Run the job
    RunJob(
      runner = runner, 
      jobname = "bascet_minhash_hist",
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
#' @inheritParams template_BascetFunction
#' @param useKMERs description
#' @return TODO
#' @export
BascetQueryFq <- function(
    bascetRoot, 
    inputName="filtered",
    outputName="kmer_counts", 
    useKMERs,
    max_reads=1000000, 
    runner=GetDefaultBascetRunner(),
    bascet_instance=GetDefaultBascetInstance()
){
  
  
  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, detect_shards_for_file(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }
  
  if(length(useKMERs)==0){
    stop("No KMERs")
  }
  
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "h5", num_shards)
  
  if(bascet_check_overwrite_output(outputFile, overwrite)) {
    #Run the job
    RunJob(
      runner = runner, 
      jobname = "bascet_query",
      cmd = c(
        shellscript_set_tempdir(bascet_instance),
        shellscript_make_bash_array("files_in",inputFiles),
        shellscript_make_bash_array("files_out",outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) helper_cancel_job_if_file_exists("${files_out[$TASK_ID]}"),
        
        shellscript_make_one_file_expander("KMERFILE", useKMERs), 
        paste(
          bascet_instance@prepend_cmd,
          bascet_instance@bin, 
          "query-fq",
          "-t $BASCET_TEMPDIR",
          "-m ", format(max_reads, scientific=FALSE),
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
#' @return KMER histogram
#' @export
BascetReadMinhashHistogram <- function(
    bascetRoot,
    inputName="minhash_hist.csv"
){
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
#' @param fname description
#' @param num_pick description
#' @param minfreq description
#' @param maxfreq description
#' @return List of KMERs
#' @export
ChooseInformativeKMERs <- function(
    kmer_hist,
    num_pick=1000, 
    minfreq=0.01, 
    maxfreq=0.10
) {
  
  ## Possibly expensive to get them all... is the total count stored somewhere? can we sample?
  kmer_hist <- kmer_hist[order(kmer_hist$cnt),]
  num_kmer <- nrow(kmer_hist)
  
  all_kmer <- kmer_hist$kmer[round(num_kmer*minfreq):round(num_kmer*maxfreq)]
  
  if(length(all_kmer)<num_pick){
    stop(paste("Cannot pick",num_pick," KMERS; only",length(all_kmer),"pass frequency criteria"))  
  }
  
  set.seed(0)
  sample(all_kmer, num_pick)
}



