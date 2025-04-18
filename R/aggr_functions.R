
################################################################################
################## Example aggregate-functions #################################
################################################################################


aggr.example <- function(
    bascetFile, 
    cellID, 
    bascet_instance
){
  
  ###### Option #1 -- Store the content in a temporary file that you have to remove once done
  tmp <- BascetReadFile(bascetFile, cellID, "out.csv", as="tempfile", bascet_instance=bascet_instance)
  my_df <- read.csv(tmp)
  file.remove(tmp)
  
  ###### Option #2 -- Get the content directly, possibly higher performance
  my_data <- BascetReadFile(bascetFile, cellID, "out.csv", as="text")   #equivalent to readLines() i.e. one big string returned (or get a list, one per line?)
  
  ###### Do something with the data and return as a data.frame
  return(data.frame(
    quality=666,
    completeness=50
  ))
}





################################################################################
################## Normal aggregate-functions ##################################
################################################################################



###############################################
#' Callback function for aggregating QUAST data.
#' To be called from BascetAggregateMap
#' 
#' @return QUAST data for each cell
#' @export
aggr.quast <- function(
    bascetFile, 
    cellID, 
    bascet_instance
){
  
  #print(cellID)
  fcont <- BascetReadFile(
    bascetFile, 
    cellID, 
    "transposed_report.tsv", 
    as="text", 
    bascet_instance=bascet_instance
  )
  if(!is.null(fcont)){
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
  } else {
    data.frame()
  }
}



###############################################
#' Callback function for aggregating QUAST data.
#' To be called from BascetAggregateMap
#' 
#' @return QUAST data for each cell
#' @export
aggr.quast_via_filesystem <- function(
    bascetFile, 
    cellID, 
    bascet_instance
){
  
  
  ### For testing
  if(FALSE){
    bascetFile <- OpenBascet(bascetRoot,"quast")
    #cellID <- "A1_B5_H8_H10"
    cellID <- "_A2_D5_H8_D12"
    bascet_instance <- bascet_instance.default
    
    foo <- BascetListFilesForCell(bascetFile,cellID)
    foo <- foo[foo$file!="cellmap.log",]
    foo <- foo[foo$cell==cellID,]
    foo
    
    #can find info here on e.g. if contigs where too short to be analyzed
    tmp <- BascetReadFile(bascetFile, cellID, "quast.log", as="tempfile", bascet_instance=bascet_instance)
    readLines(tmp)
  }
  
  #print(cellID)
  tmp <- BascetReadFile(
    bascetFile, 
    cellID, 
    "transposed_report.tsv", 
    as="tempfile", 
    bascet_instance=bascet_instance
  )
  if(!is.null(tmp)){
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
  } else {
    data.frame()
  }
}









###############################################
#' Callback function for aggregating min-hashes for each cell
#' To be called from BascetAggregateMap
#' 
#' @return Minhash data (minhash.txt) for each cell
#' @export
aggr.minhash <- function(
    bascetFile, 
    cellID, 
    bascet_instance
){
  dat <- BascetReadFile(bascetFile, cellID, "minhash.txt", as="text", bascet_instance=bascet_instance)
#  dat <- readLines(tmp)
#  file.remove(tmp)
  
  set_kmer <- stringr::str_split_i(dat,"\t",1)
  set_kmer
}





###############################################
#' Callback function for aggregating min-hashes for each cell
#' To be called from BascetAggregateMap
#' 
#' @return Minhash data (minhash.txt) for each cell
#' @export
aggr.minhash_via_fs <- function(
    bascetFile, 
    cellID, 
    bascet_instance
){
  tmp <- BascetReadFile(bascetFile, cellID, "minhash.txt", as="tempfile", bascet_instance=bascet_instance)
  dat <- readLines(tmp)
  file.remove(tmp)
  
  set_kmer <- stringr::str_split_i(dat,"\t",1)
  set_kmer
}



################################################################################
################## Simplified calls to aggregate ###############################
################################################################################


###############################################
#' Aggregate frequency of minhashes across cells
#' 
#' @inheritParams template_BascetFunction
#' @param inputName description
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



################################################################################
################## Non-standard aggregate-functions ############################
################################################################################




###############################################
#' Callback function for just getting raw file contents
#' To be called from BascetAggregateMap
#' 
#' @param fname Name of output filename to get for each cell
#' @return QUAST data for each cell
#' @export
aggr.rawtext <- function(fname){
  function(bascetFile, cellID, bascet_instance){
    rawtext <- BascetReadFile(bascetFile, cellID, fname, as="text")  
    data.frame(
      rawtext=rawtext
    )
  }
}

