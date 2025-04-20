


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
      
      ### Abort early if needed    
      if(!overwrite) helper_cancel_job_if_file_exists("${files_out[$TASK_ID]}"),

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
      jobname = paste0("Z_map_",withfunction),
      cmd = cmd,
      arraysize = num_shards
    )      
  } else {
    new_no_job()
  }
}




###############################################
#' Convenience function; alternative is to somehow implement as.data.frame
#' @return TODO
#' @export
MapListAsDataFrame <- function(mylist){
  
  if(FALSE){
    mylist <- list()
    mylist[["a"]] <- data.frame(x=6)
    mylist[["b"]] <- data.frame(x=3)
#    mylist[["c"]] <- data.frame()
    mylist[["c"]] <- NULL
  }
  
  #Make table
  out <- do.call(rbind, mylist)

  #Some functions return multiple, or no, lines. figure out count
  num_entry <- sapply(mylist, nrow)
  final_row_name <- rep(-1, nrow(out))
  cur_row <- 1
  for(i in 1:length(num_entry)){
    final_row_name[cur_row:(cur_row-1+num_entry[i])] <- i
    cur_row <- cur_row+num_entry[i]
  }
  
  #Set row names based on index
  rownames(out) <- names(mylist)[final_row_name]    #names(mylist)[!is.null(mylist)]
  out
}






###############################################
#' Convenience function; alternative is to somehow implement as.data.frame
#' 
#' This one puts cellID as a new column
#' 
#' make this the new MapListAsDataFrame ?
#' 
#' @return TODO
#' @export
MapCellMultiListAsDataFrame <- function(mylist){
  
  if(FALSE){
    mylist <- list()
    mylist[["a"]] <- data.frame(x=6)
    mylist[["b"]] <- data.frame(x=3)
    mylist[["c"]] <- NULL
  }
  
  newlist <- list()
  for(i in names(mylist)){
    one_entry <- mylist[[i]]
    if(nrow(one_entry)>0){
      one_entry$cellID <- i
      newlist[[i]] <- one_entry
    }
  }
  
  out <- do.call(rbind, newlist)
  
  out
}



###############################################
#' Count entries in long format data frame and return as a sparse matrix
#' 
#' @param dat A data.frame
#' @param rowname Column to use as matrix row name
#' @param colname Column to use as matrix column name
#' @return A dgCMatrix sparse matrix
#' @export
CountDataFrameToSparseMatrix <- function(dat, rowname, colname) {
  red_dat <- data.frame(
    col=dat[,colname],
    row=dat[,rowname]
  )
  red_dat <- sqldf::sqldf("select col, row, count(*) as cnt from red_dat group by col, row") 
  
  red_dat$col <- factor(red_dat$col)
  red_dat$row <- factor(red_dat$row)
  mat <- Matrix::sparseMatrix(  
    i=as.integer(red_dat$row),
    j=as.integer(red_dat$col), 
    x=red_dat$cnt
  )
  colnames(mat) <- levels(red_dat$col)
  rownames(mat) <- levels(red_dat$row)
  mat
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
    include_cells=NULL,
    showProgress=TRUE,
    verbose=FALSE,
    bascet_instance=GetDefaultBascetInstance()
){
  
  #Get file coordinates of all objects in zip file
  cellname_coord <- BascetCellNames(bascetRoot, bascetName)
  
  #Open the file, prep for reading
  if(verbose){
    print("Creating extract streamer session")
  }
  bascetFile <- OpenBascet(bascetRoot, bascetName)
  if(verbose){
    print("Extract streamer session ok")
  }
  pbar <- progress::progress_bar$new(total = length(bascetFile@cellmeta$cell))
  if(showProgress){
    pbar$tick(0)
  }
  
  #Loop over all cells by default
  if(is.null(include_cells)){
    include_cells <- bascetFile@cellmeta$cell
  }
  
  #Loop over all files in the bascet
  output <- list()
  for(cellname in include_cells){
    output[[cellname]] <- aggrFunction(bascetFile, cellname, bascet_instance)
    if(showProgress){
      pbar$tick()
    }
  }
  
  #End associated bascet session
  CloseBascet(bascetFile)
  
  output
}


















################################################################################
########### Testing ############################################################
################################################################################


if(FALSE){
  bascetRoot <- "/home/mahogny/jupyter/bascet/zorn/try_unzip"
  bascetName <- "quast"
  #aggrFunction <- aggr.quast
  BascetAggregateMap("/home/mahogny/jupyter/bascet/zorn/try_unzip","quast",aggr.quast) 
  
  #tmp <- BascetReadMapFile(bascetFile, cellID, "out.csv", as="tempfile")

  
  
  
  fquast <- OpenBascet(bascetRoot, "quast")
  allf <- BascetListFilesForCell(fquast, "D1_B4_F7_B12")               ## todo for all cells, and one cell!
  allf
  
  allf[allf$file!="cellmap.log",]

  tfile <- BascetReadFile(fskesa, "E2_B4_E9_E11", "transposed_report.tsv", as="tempfile")
  readLines(tfile)
  
  
}







