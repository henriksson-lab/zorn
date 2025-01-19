
### Where all data is stored
bascetRoot = "testdata"  


### Run 
BascetMapCell(
  bascetRoot,
  withfunction = "./my_map_script.sh",
  inputName = "filtered",
  outputName = "my_output",
  runner=inst
)





######### Callback function for aggregating minhash
aggr.myfunc <- function(bascetFile, cellID){
  tmp <- BascetReadFile(bascetFile, cellID, "minhash.txt", as="tempfile")
  dat <- readLines(tmp)
  file.remove(tmp)

  set_kmer <- stringr::str_split_i(dat,"\t",1)
  set_kmer
}





my_aggr <- BascetAggregateMap(bascetRoot,"my_output",aggr.myfunc) 
