



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







######### Callback function for aggregating QUAST data
aggr.quast <- function(bascetFile, cellID){
  ###### Option #1
  tmp <- BascetReadFile(bascetFile, cellID, "transposed_report.tsv", as="tempfile")
  dat <- readLines(tmp)
  file.remove(tmp)
  
  
  return(data.frame(
    quality=666,
    completeness=50
  ))
}




if(FALSE){
  
  
  fquast <- OpenBascet(bascetRoot, "quast")
  allf <- BascetListFilesForCell(fquast, "D1_B4_F7_B12")               ## todo for all cells, and one cell!
  allf
  
  allf[allf$file!="cellmap.log",]

  tfile <- BascetReadFile(fskesa, "E2_B4_E9_E11", "transposed_report.tsv", as="tempfile")
  readLines(tfile)
  
  
}

