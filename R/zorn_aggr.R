



######### Example callback function for aggregating data
aggr.example <- function(bascetFile, cellID){
  ###### Option #1
  tmp <- BascetReadFile(bascetFile, cellID, "out.csv", as="tempfile")
  my_df <- read.csv(tmp)
  file.remove(tmp)
  
  ###### Option #2
  my_data <- BascetReadFile(bascetFile, cellID, "out.csv", as="text")   #equivalent to readLines() i.e. one big string returned (or get a list, one per line?)
  
  
  return(data.frame(
    quality=666,
    completeness=50
  ))
}




################################################################################
################## Aggregate-functions #########################################
################################################################################



######### Callback function for aggregating QUAST data
#' @return TODO
#' @export
aggr.quast <- function(bascetFile, cellID){
  ###### Option #1
  tmp <- BascetReadFile(bascetFile, cellID, "transposed_report.tsv", as="tempfile")
  dat <- readLines(tmp)
  file.remove(tmp)
  
  
  #TODO got lost, need to implemented again
  
  return(data.frame(
    quality=666,
    completeness=50
  ))
}









######### Callback function for aggregating minhash
#' @return TODO
#' @export
aggr.minhash <- function(bascetFile, cellID){
  tmp <- BascetReadFile(bascetFile, cellID, "minhash.txt", as="tempfile")
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
    inputName="minhash"
) {
  minhash_aggr <- BascetAggregateMap(
    bascetRoot,
    inputName,
    aggr.minhash
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

