





######### Example callback function for aggregating data
aggr.quast <- function(bascetFile, cellID){
  ###### Option #1
  tmp <- BascetReadMapFile(bascetFile, cellID, "out.csv", as="tempfile")  #TODO: could support on the fly conversion, FASTQ => FASTA, as needed
  my_df <- read.csv(tmp)
  file.remove(tmp)
  
  ###### Option #2
  my_data <- BascetReadMapFile(bascetFile, cellID, "out.csv", as="text")   #equivalent to readLines() i.e. one big string returned (or get a list, one per line?)
  
  
  return(data.frame(
    quality=666,
    completeness=50
  ))
}

