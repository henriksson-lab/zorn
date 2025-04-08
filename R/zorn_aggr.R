


################################################################################
################ The map system ################################################
################################################################################


#had a copy in zorn.R; check there if there are issues

###############################################
#' Call a MAP function for all cells
#' 
#'     TODO override args docs
#' 
#' @param args List of arguments (key,value) to provide to the script
#' 
#' @export
BascetMapCell <- function(
    bascetRoot, 
    withfunction, 
    inputName, 
    outputName, 
    args=list(),
    overwrite=TRUE,
    runner,
    bascet_instance=bascet_instance.default
){
  
  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, detect_shards_for_file(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }
  
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "zip", num_shards)
  
  if(bascet_check_overwrite_output(outputFiles, overwrite)) {
    #Build the command - custom arguments
    cmd <- c()
    for(key in names(args)){
      
      #Escape value
      val <- args[[key]]
      val <- stringr::str_replace_all(val,stringr::fixed("\""),"\\\"")
      
      #Add argument
      cmd <- c(
        cmd,
        paste0("export ",key,"=\"",val,"\"")
      )
    }
    
    #Build the command - the rest
    cmd <- c(
      cmd,
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
    )
    
    #Run the job
    RunJob(
      runner = runner, 
      jobname = paste0("bascet_map_",withfunction),
      cmd = cmd,
      arraysize = num_shards
    )      
  }
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
 
  
  
  
  fquast <- OpenBascet(bascetRoot, "quast")
  allf <- BascetListFilesForCell(fquast, "D1_B4_F7_B12")               ## todo for all cells, and one cell!
  allf
  
  allf[allf$file!="cellmap.log",]

  tfile <- BascetReadFile(fskesa, "E2_B4_E9_E11", "transposed_report.tsv", as="tempfile")
  readLines(tfile)
  
  
}

