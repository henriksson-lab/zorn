



###############################################
#' Select kmers that appear useful for clustering
#' 
#' 
#' CHANGED: it just sums up kmers
#' 
#' @return TODO
#' @export
BascetFeaturise <- function( ########### need a better name; KMC something?
    bascetRoot, 
    inputName="kmc", 
    outputName="sumkmers", 
    includeCells=NULL,
    runner,
    bascet_instance=bascet_instance.default){
  
  
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
  
  #Run the job
  RunJob(
    runner = runner, 
    jobname = "bascet_query",
    cmd = c(
      shellscript_set_tempdir(bascet_instance),
      if(produce_cell_list) shellscript_make_files_expander("CELLFILE", list_cell_for_shard),
      shellscript_make_bash_array("files_in",inputFiles),
      shellscript_make_bash_array("files_out",outputFiles),
      paste(
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
  
}






###############################################
#' Build count table from kmer table and selected kmers
#' @return TODO
#' @export
BascetQuery <- function(
    bascetRoot, 
    inputName="kmc",
    outputName="kmer", 
    useKMERs,
    runner,
    bascet_instance=bascet_instance.default){
  
  
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
  
  #Run the job
  RunJob(
    runner = runner, 
    jobname = "bascet_query",
    cmd = c(
      shellscript_set_tempdir(bascet_instance),
      shellscript_make_bash_array("files_in",inputFiles),
      shellscript_make_bash_array("files_out",outputFiles),
      shellscript_make_one_file_expander("KMERFILE", useKMERs), 
      paste(
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
  
}







################################################################################
########################## Count matrix tools ##################################
################################################################################


###############################################
#' Read a count matrix as produced by Bascet (hdf5 format)
#' @return TODO
#' @export
ReadBascetCountMatrix <- function(fname){
  print("Loading HDF5 file")
  h5f <- rhdf5::H5Fopen(fname)
  indices <- h5f$X$indices+1
  indptr <-  h5f$X$indptr
  dat <- h5f$X$data
  shape <- h5f$X$shape

  print(paste0("Assembling matrix, size: ", shape[1],"x",shape[2]))
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
############################# KMER tools #######################################
################################################################################


###############################################
#' Generate a histogram from a KMC3 database
#' @return TODO
#' @export
KmcGetHistogram <- function(fname, plot=TRUE){
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
#' @return TODO
#' @export
KmcGetKmers <- function(fname, mincount=NULL, maxcount=NULL) {
  
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
#' @return TODO
#' @export
KmcChooseKmerFeatures <- function(fname, num_pick=1000, minfreq=0.01, maxfreq=0.10) {
  
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
