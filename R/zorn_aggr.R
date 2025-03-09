


################################################################################
################ The map system ################################################
################################################################################



###############################################
#' Call a function for all cells
#' @return TODO
#' @export
BascetMapCell <- function(
    bascetRoot, 
    withfunction, 
    inputName, 
    outputName, 
    runner, 
    bascet_instance=bascet_instance.default){
  
  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, detect_shards_for_file(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }
  
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "zip", num_shards)
  
  #Run the job
  RunJob(
    runner = runner, 
    jobname = paste0("bascet_map_",withfunction),
    cmd = c(
      shellscript_set_tempdir(bascet_instance),
      shellscript_make_bash_array("files_in",inputFiles),
      shellscript_make_bash_array("files_out",outputFiles),
      paste(
        bascet_instance@prepend_cmd,
        bascet_instance@bin, 
        "mapcell",
        "-t $BASCET_TEMPDIR",
        "-i ${files_in[$TASK_ID]}",
        "-o ${files_out[$TASK_ID]}",
        "-s", withfunction)
    ),
    arraysize = num_shards
  )  
}




###############################################
#' Convenience function; alternative is to somehow implement as.data.frame
#' @return TODO
#' @export
MapListAsDataFrame <- function(mylist){
  out <- do.call(rbind, mylist)
  rownames(out) <- names(mylist)
  out
}




###############################################
#' Aggregate data from previous Map call
#' 
#' todo: allow multi-cpu support? parallel library
#
#' todo note, this is effectively a pure-R map function. different name?
#' @return TODO
#' @export
BascetAggregateMap <- function(
    bascetRoot, 
    bascetName, 
    aggrFunction,
    showProgress=TRUE,
    bascet_instance
){
  
  #Get file coordinates of all objects in zip file
  cellname_coord <- BascetCellNames(bascetRoot, bascetName)
  
  #Open the file, prep for reading
  bascetFile <- OpenBascet(bascetRoot, bascetName)
  
  pbar <- progress::progress_bar$new(total = length(bascetFile@cellmeta$cell))
  if(showProgress){
    pbar$tick(0)
  }
  
  #Loop over all files in the bascet
  output <- list()
  for(cellname in bascetFile@cellmeta$cell){
    output[[cellname]] <- aggrFunction(bascetFile, cellname, bascet_instance)
    if(showProgress){
      pbar$tick()
    }
  }
  output
}



if(FALSE){
  bascetRoot <- "/home/mahogny/jupyter/bascet/zorn/try_unzip"
  bascetName <- "quast"
  aggrFunction <- aggr.quast
  BascetAggregateMap("/home/mahogny/jupyter/bascet/zorn/try_unzip","quast",aggr.quast) 
  
  #tmp <- BascetReadMapFile(bascetFile, cellID, "out.csv", as="tempfile")
  
}






################################################################################
################## Aggregate-functions #########################################
################################################################################






######### Example callback function for aggregating data
aggr.example <- function(bascetFile, cellID, bascet_instance){
  ###### Option #1
  tmp <- BascetReadFile(bascetFile, cellID, "out.csv", as="tempfile", bascet_instance=bascet_instance)
  my_df <- read.csv(tmp)
  file.remove(tmp)
  
  ###### Option #2
  my_data <- BascetReadFile(bascetFile, cellID, "out.csv", as="text")   #equivalent to readLines() i.e. one big string returned (or get a list, one per line?)
  
  
  return(data.frame(
    quality=666,
    completeness=50
  ))
}







######### Callback function for aggregating QUAST data
#' @return TODO
#' @export
aggr.quast <- function(bascetFile, cellID, bascet_instance){
  
  if(FALSE){
    bascetFile <- OpenBascet(bascetRoot,"quast")
    #"quast.1.zip"
    #    bascetFile()
    cellID <- "A1_B5_H8_H10"
    bascet_instance <- bascet_inst
    
    foo <- BascetListFilesForCell(bascetFile,cellID)
    foo <- foo[foo$file!="cellmap.log",]
    foo
  }
  
  
  if(FALSE){
    #can find info here on e.g. if contigs where too short to be analyzed
    tmp <- BascetReadFile(bascetFile, cellID, "quast.log", as="tempfile", bascet_instance=bascet_instance)
    readLines(tmp)
  }
  
  if(FALSE){
    #General output of the script
    tmp <- BascetReadFile(bascetFile, cellID, "cellmap.log", as="tempfile", bascet_instance=bascet_instance)
    readLines(tmp)
  }  
  
  
  
  #tmp <- BascetReadFile(bascetFile, cellID, "transposed_report2.tsv", as="tempfile", bascet_instance=bascet_instance) #how to detect error?
  tmp <- BascetReadFile(bascetFile, cellID, "transposed_report.tsv", as="tempfile", bascet_instance=bascet_instance)
  fcont <- readLines(tmp, n=2)
  #dat <- read.table(tmp)
  file.remove(tmp)
  
  dat <- data.frame(
    row.names=stringr::str_split(fcont[1],"\t")[[1]],
    value=stringr::str_split(fcont[2],"\t")[[1]]
  )
  dat <- dat[-1,,drop=FALSE]
  
  rownames(dat) <- stringr::str_replace_all(rownames(dat), stringr::fixed("#"),"Number of")
  
  #Arrange in the right format
  dat <- t(dat)
  
  #TODO set data types to double whenever possible
  
  dat
}








######### Callback function for aggregating minhash
#' @return TODO
#' @export
aggr.minhash <- function(bascetFile, cellID, bascet_instance){
  tmp <- BascetReadFile(bascetFile, cellID, "minhash.txt", as="tempfile", bascet_instance=bascet_instance)
  dat <- readLines(tmp)
  file.remove(tmp)
  
  set_kmer <- stringr::str_split_i(dat,"\t",1)
  set_kmer
}



################################################################################
################## Simplified calls to aggregate ###############################
################################################################################



#' Aggregate frequency of minhashes across cells
#' @return TODO
#' @export
AggregateMinhashes <- function(
    bascetRoot,
    inputName="minhash",
    bascet_instance
) {
  minhash_aggr <- BascetAggregateMap(
    bascetRoot,
    inputName,
    aggr.minhash,
    bascet_instance=bascet_instance
  )
  
  #Put KMERs into a data frame
  all_kmer <- do.call(c, minhash_aggr)
  all_kmer <- as.data.frame(table(all_kmer))
  colnames(all_kmer) <- c("kmer","freq")
  
  #Order KMERs
  all_kmer <- all_kmer[order(all_kmer$freq, decreasing = TRUE),]
  all_kmer$index <- 1:nrow(all_kmer)
  
  all_kmer
}








if(FALSE){
  
  
  fquast <- OpenBascet(bascetRoot, "quast")
  allf <- BascetListFilesForCell(fquast, "D1_B4_F7_B12")               ## todo for all cells, and one cell!
  allf
  
  allf[allf$file!="cellmap.log",]

  tfile <- BascetReadFile(fskesa, "E2_B4_E9_E11", "transposed_report.tsv", as="tempfile")
  readLines(tfile)
  
  
}

