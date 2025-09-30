################################################################################
################ The map system ################################################
################################################################################

###############################################
#' Call a MAP function for all cells
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param withfunction Function to apply
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param args List of arguments (key,value) to provide to the script
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetMapCell <- function(
    bascetRoot, 
    withfunction, 
    inputName, 
    outputName, 
    args=list(),
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  #check arguments
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.function(withfunction))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.list(args))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  
  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, detectShardsForFile(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }
  
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "zip", num_shards)
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
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
      #shellscript_set_tempdir(bascetInstance),
      shellscriptMakeBashArray("files_in",inputFiles),
      shellscriptMakeBashArray("files_out",outputFiles),
      
      ### Abort early if needed    
      if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),

      paste(
        bascetInstance@prependCmd,
        bascetInstance@bin, 
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
      bascetInstance = bascetInstance,
      cmd = cmd,
      arraysize = num_shards
    )      
  } else {
    new_no_job()
  }
}




###############################################
#' Convenience function; alternative is to somehow implement as.data.frame.
#' 
#' @param mylist TODO
#' 
#' @return A data.frame
#' @export
MapListAsDataFrame <- function(
    mylist
){
  #check arguments
  stopifnot(is.list(mylist))
  
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
  rownames(out) <- names(mylist)[final_row_name]
  #names(mylist)[!is.null(mylist)]
  out
}






###############################################
#' Convenience function; alternative is to somehow implement as.data.frame
#' 
#' This one puts cellID as a new column
#' 
#' FUTURE possible make this the new MapListAsDataFrame
#' 
#' @return TODO
#' @export
MapCellMultiListAsDataFrame <- function(mylist){
  #check arguments
  stopifnot(is.list(mylist))
  
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
  
  do.call(rbind, newlist)
}



###############################################
#' Count entries in long format data frame and return as a sparse matrix
#' 
#' @param dat A data.frame
#' @param rowname Column to use as matrix row name
#' @param colname Column to use as matrix column name
#' 
#' @return A dgCMatrix sparse matrix
#' @export
CountDataFrameToSparseMatrix <- function(
    dat, 
    rowname, 
    colname
) {
  #check arguments
  stopifnot(is.data.frame(dat))
  stopifnot(is.character(colname))
  stopifnot(is.character(rowname))

  #Extract relevant subset of dataframe  
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
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param aggrFunction Function to use for extracting simplified data for cell
#' @param includeCells Cells to aggregate
#' @param showProgress Show progress bar
#' @param verbose Show debug output
#' 
#' @return Aggregated data
#' @export
BascetAggregateMap <- function(
    bascetRoot, 
    inputName, 
    aggrFunction,
    includeCells=NULL,
    showProgress=TRUE,
    verbose=FALSE,
    bascetInstance=GetDefaultBascetInstance()
){
  #check arguments
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.function(aggrFunction))
  stopifnot(is.valid.listcells(includeCells))
  stopifnot(is.logical(showProgress))
  stopifnot(is.logical(verbose))
  stopifnot(is.bascet.instance(bascetInstance))
  
  
  #Get file coordinates of all objects in zip file
  cellname_coord <- BascetCellNames(bascetRoot, inputName, bascetInstance)  ############## todo: avoid opening streamer twice
  
  #Open the file, prep for reading
  if(verbose){
    print("Creating extract streamer session")
  }
  bascetFile <- OpenBascet(bascetRoot, inputName)
  if(verbose){
    print("Extract streamer session ok")
  }
  pbar <- progress::progress_bar$new(total = length(bascetFile@cellmeta$cell))
  if(showProgress){
    pbar$tick(0)
  }
  
  #Loop over all cells by default
  if(is.null(includeCells)){
    includeCells <- bascetFile@cellmeta$cell
  }
  
  #Loop over all files in the bascet
  output <- list()
  for(cellname in includeCells){
    output[[cellname]] <- aggrFunction(bascetFile, cellname, bascetInstance)
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
