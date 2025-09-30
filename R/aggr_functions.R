
################################################################################
################## Example aggregate-functions #################################
################################################################################


aggr.example <- function(
    bascetFile, 
    cellID, 
    bascetInstance
){
  
  if(TRUE){
    ###### Option #1 -- Store the content in a temporary file that you have to remove once done
    tmp <- BascetReadFile(bascetFile, cellID, "out.csv", as="tempfile", bascetInstance=bascetInstance)
    my_df <- read.csv(tmp)
    file.remove(tmp)
    
  }

  if(TRUE){
    ###### Option #2 -- Get the content directly, possibly higher performance
    my_data <- BascetReadFile(bascetFile, cellID, "out.csv", as="text")   #equivalent to readLines() i.e. one big string returned (or get a list, one per line?)
    
    # If you do option #2, you can do this to parse the lines using common table readers
    zz <- textConnection(my_data)
    dat <- read.delim(zz)
    close(zz)    
  }  

  ###### Do something with the data and return ideally as a one-line data.frame (multiple lines if multiple outputs).
  ###### Note #1: there is no notion of cellIDs here. this is handled on a higher abstraction level
  ###### Note #2: you can return the data in any shape you want. just be prepared to handle it after aggregation!
  return(data.frame(
    quality=666,
    completeness=50
  ))
}


################################################################################
################## Core aggregate-functions ####################################
################################################################################



###############################################
#' Callback function for aggregating min-hashes for each cell.
#' To be called from BascetAggregateMap
#' 
#' @param bascetFile Bascet file handle
#' @param cellID ID of cell to process
#' @param bascetInstance A Bascet instance
#' 
#' @return Minhash data (minhash.txt) for each cell
#' @export
aggr.minhash <- function(
    bascetFile, 
    cellID, 
    bascetInstance
){
  dat <- BascetReadFile(bascetFile, cellID, "minhash.txt", as="text", bascetInstance=bascetInstance)

  set_kmer <- stringr::str_split_i(dat,"\t",1)
  set_kmer
}




################################################################################
################## Simplified calls to aggregate ###############################
################################################################################


###############################################
#' Aggregate frequency of minhashes across cells
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard (Container of minhashes)
#' @param bascetInstance A Bascet instance
#' 
#' @return Data.frame of KMERs and frequencies
#' @export
AggregateMinhashes <- function(
    bascetRoot,
    inputName="minhash",
    bascetInstance
) {
  minhash_aggr <- BascetAggregateMap(
    bascetRoot,
    inputName,
    aggr.minhash,
    bascetInstance=bascetInstance
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
#' 
#' @return QUAST data for each cell
#' @export
aggr.rawtext <- function(
    fname
){
  function(bascetFile, cellID, bascetInstance){
    rawtext <- BascetReadFile(bascetFile, cellID, fname, as="text")  
    
    data.frame(
      rawtext=rawtext
    )
  }
}

