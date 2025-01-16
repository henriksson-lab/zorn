
if(FALSE){
  #to run roxygen --- run in clean workspace!
  devtools::document()
}

################################################################################
################ Settings for the underlying Bascet installation ###############
################################################################################

#' @export
setClass("BascetInstance", slots=list(
  bin="character",
  tempdir="character"
)
) 


###############################################
#' Create a new bascet instance
#' @return TODO
#' @export
BascetInstance <- function(bin, tempdir){
  if(!is.null(tempdir) & !file.exists(tempdir)){
    stop(sprintf("temp directory %s does not exist", tempdir))
  }
  new(
    "BascetInstance",
    bin=bin,
    tempdir=tempdir
  )
}



###############################################
#' The default Bascet installation settings
#' @return TODO
#' @export
bascet_instance.default <- BascetInstance(
  bin="/home/mahogny/jupyter/bascet/target/debug/robert",
  tempdir="/data/henlab/bascet_temp"
)

#bascet_instance.default <- BascetInstance(bin="bascet")




###############################################
#' Get a temp directory to use; need to be created
#' @return TODO
#' @export
GetBascetTempDir <- function(bascet_instance){
  if(is.null(bascet_instance@tempdir)){
    tempfile()
  } else {
    tempfile(tmpdir=file.path(bascet_instance@tempdir))
  }
}





################################################################################
################ Bascet command line tool: debarcoding raw fastq ###############
################################################################################



###############################################
#' Detect metadata for raw input FASTQ files
#' 
#' @param rawRoot Path to folder with FASTQ files
#' @return A data frame with metadata for the raw input files
#' @examples
#' DetectRawFileMeta("/path/to/raw_fastq", verbose = TRUE)
DetectRawFileMeta <- function(rawRoot, verbose=FALSE){
  #rawRoot <- "/husky/fromsequencer/241206_novaseq_wgs3/raw"
  allfiles <- list.files(rawRoot)
  
  if(verbose){
    print("Found the following files in the directory")
    print(allfiles)
  }
  
  allfiles <- allfiles[
    stringr::str_ends(allfiles,".fq.gz") | 
      stringr::str_ends(allfiles, ".fastq.gz") | 
      stringr::str_ends(allfiles, ".fq") | 
      stringr::str_ends(allfiles,".fastq")]
  
  if(verbose){
    print("Found the following raw read-like files")
    print(allfiles)
  }
  
  r1_files <- allfiles[stringr::str_detect(allfiles,stringr::fixed("_R1"))]
  
  if(verbose){
    print("Found the following R1-like files")
    print(r1_files)
  }
  
  r2_corresponding <- stringr::str_replace(r1_files,stringr::fixed("R1"),"R2")
  
  if(!all(r2_corresponding %in% allfiles)){
    stop("Not all R1 files appear to have a corresponding R2 file")
  }
  
  ### Gather first set of metadata
  meta <- data.frame(
    r1=r1_files,
    r2=r2_corresponding,
    dir=file.path(rawRoot) #this standarizes the trailing /
  )
  meta$prefix <- ""
  
  #Guess prefixes
  meta$possible_prefix <- stringr::str_split_i(meta$r1, "_S[0123456789]",1)
  unique_prefix <- unique(meta$possible_prefix)
  if(verbose){
    print("Detected possible prefixes")
    print(unique_prefix)
  }
  
  #If there is more than one prefix, then we have to add them. otherwise just keep it simple
  if(length(unique_prefix)>1){
    meta$prefix <- paste0(meta$possible_prefix,"_")
  } else {
    if(verbose){
      print("Detected one one possible prefix, so not adding it")
    }
  }
  
  #Sanitize prefixes. Some characters will break BAM tags etc
  meta$prefix <- stringr::str_remove_all(meta$prefix, " ")
  meta$prefix <- stringr::str_remove_all(meta$prefix, "/")
  meta$prefix <- stringr::str_remove_all(meta$prefix, "\"")
  
  #Return metadata
  meta[,c("prefix","r1","r2","dir")]
}

#rawmeta <- DetectRawFileMeta("/husky/fromsequencer/241206_novaseq_wgs3/raw")


###############################################
#' Generate BAM with barcodes from input raw FASTQ
#' @return TODO
#' @export
BascetGetRawAtrandiWGS <- function(
    bascetRoot, 
    rawmeta, 
    outputName="debarcoded", 
    outputNameIncomplete="incomplete_reads", 
    chemistry="atrandi_wgs",  #or atrandi_rnaseq
    runner, 
    bascet_instance=bascet_instance.default){

  #One output per input pair of reads  
  num_shards <- nrow(rawmeta)
  if(num_shards==0){
    stop("No input files")
  }
  outputFilesComplete <- make_output_shard_names(bascetRoot, outputName, "tirp.gz", num_shards)
  outputFilesIncomplete <- make_output_shard_names(bascetRoot, outputNameIncomplete, "tirp.gz", num_shards)
  

  RunJob(
    runner = runner, 
    jobname = "bascet_getraw",
    cmd = c(
      shellscript_set_tempdir(bascet_instance),
      shellscript_make_bash_array("files_r1",file.path(rawmeta$dir, rawmeta$r1)),
      shellscript_make_bash_array("files_r2",file.path(rawmeta$dir, rawmeta$r2)),
      shellscript_make_bash_array("files_out",outputFilesComplete),
      shellscript_make_bash_array("files_out_incomplete",outputFilesIncomplete),
      paste(
        bascet_instance@bin, 
        "getraw",
        "-t $BASCET_TEMPDIR",
        "--chemistry",chemistry,  
        "--r1 ${files_r1[$TASK_ID]}",  
        "--r2 ${files_r2[$TASK_ID]}",
        "--out-complete   ${files_out[$TASK_ID]}",                 #Each job produces a single output
        "--out-incomplete ${files_out_incomplete[$TASK_ID]}"                 #Each job produces a single output
      )
    ),
    arraysize = nrow(rawmeta)
  )
}




################################################################################
################ Bascet command line tools: RNA-seq ############################
################################################################################


###############################################
#' Aligned debarcoded BAMs
#'  #sort or not here?
#' @return TODO
#' @export
BascetAlign <- function(bascetRoot, genomeReference, inputName="debarcoded", outputName="unsorted_aligned", runner, bascet_instance=bascet_instance.default){ 
  
  #Figure out input and output file names  
  inputFiles <- detect_shards_for_file(bascetRoot, inputName, "bam")
  num_shards <- length(inputFiles)
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "bam", num_shards)
  
  #Run the job
  RunJob(
    runner = runner, 
    jobname = "bascet_align",
    withdata = c(
      "-t $BASCET_TEMPDIR",
      shellscript_set_tempdir(bascet_instance),
      shellscript_make_bash_array("files_in",inputFiles),
      shellscript_make_bash_array("files_out",outputFiles)
    ),
    cmd = paste(bascet_instance@bin, "align --in $files_in[$TASK_ID] --o $files_out[$TASK_ID] --reference ", genomeReference),
    arraysize = num_shards
  )
  
  
}


################################################################################
################ Bascet command line tools: WGS ################################
################################################################################




###############################################
#' Take debarcoded reads and split them into suitable numbers of shards
#' @return TODO
#' @export
BascetShardify <- function(
    bascetRoot, 
    inputName="debarcoded", 
    includeCells=NULL,
    num_output_shards=1,
    outputName="filtered", 
    runner, 
    bascet_instance=bascet_instance.default){

  input_shards <- detect_shards_for_file(bascetRoot, inputName)
  if(length(input_shards)==0){
    stop("Found no input files")
  }

  #Include all cells if nothing else provided
  if(is.null(includeCells)){
    includeCells <- BascetCellNames(bascetRoot, inputName)
    includeCells <- unique(includeCells$cell) #when shardifying, we expect cells to appear more than once  -- could warn for other commands!
    print(paste("Including all the",length(includeCells)), "cells")
  }

  #Figure out which cell goes into which shard
  list_cell_for_shard <- shellscript_split_arr_into_list_randomly(includeCells, num_output_shards)
  
  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, input_shards)
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "tirp.gz", num_output_shards)

  #Produce the script and run the job
  RunJob(
    runner = runner, 
    jobname = "bascet_shardify",
    cmd = c(
      shellscript_set_tempdir(bascet_instance),
      shellscript_make_files_expander("CELLFILE", list_cell_for_shard),
      shellscript_make_bash_array("files_out", outputFiles),
      paste(
        bascet_instance@bin,
        "shardify", 
        "-t $BASCET_TEMPDIR",
        "-i",shellscript_make_commalist(inputFiles), #Need to give all input files for each job
        "-o ${files_out[$TASK_ID]}",                 #Each job produces a single output
        "--cells $CELLFILE"                          #Each job takes its own list of cells
      ),  
      "rm $CELLFILE"
      ), 
    arraysize = num_output_shards
  )
  
}




###############################################
#' Assemble the genomes
#' @return TODO
#' @export
BascetAssemble <- function(bascetRoot, inputName="rawreads", outputName="assembled", runner, bascet_instance=bascet_instance.default){
  stop("this is done using MapCell; but we could produce a wrapper for it")
}




###############################################
#' Transform: subset, convert, merge, divide
#' @return TODO
#' @export
BascetMapTransform <- function(
    bascetRoot, 
    inputName, 
    outputName,
    num_divide=1,
    num_merge=1,
    out_format="tirp.gz", ### not really!
    includeCells=NULL,
    runner, 
    bascet_instance=bascet_instance.default){
  
  
  is.integer.overequal1 <- function(x){
    return(x>=1 & x==as.integer(x))
  }
  
  #Verify input
  if(!(num_divide==1 || num_merge==1)){
    stop("Either divide or merge must be set to 1")
  }
  if(!is.integer.overequal1(num_divide) || !is.integer.overequal1(num_merge)){
    stop("Number of divisions and merges must be an integer >= 1")
  }
  
  #TODO check outformat

  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, detect_shards_for_file(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }
  
  outputFiles <- make_output_shard_names(bascetRoot, outputName, out_format, num_shards)

  #If cell list is provided, produce a file for input (not all transform calls can handle this, so optional)
  produce_cell_list <- !is.null(includeCells)
  if(produce_cell_list) {
    #Figure out which cell goes into which shard   TODO, using the same cell list for all shards
    list_cell_for_shard <- list()
    for(i in 1:length(inputFiles)){
      list_cell_for_shard[[i]] <- includeCells
    }
  }

  #Run the job
  RunJob(
    runner = runner, 
    jobname = paste0("bascet_transform"),
    cmd = c(
      shellscript_set_tempdir(bascet_instance),
      if(produce_cell_list) shellscript_make_files_expander("CELLFILE", list_cell_for_shard),
      shellscript_make_bash_array("files_in",inputFiles),
      shellscript_make_bash_array("files_out",outputFiles),
      paste(
        bascet_instance@bin, 
        "transform",
        if(produce_cell_list) "--cells $CELLFILE",
        #"-t $BASCET_TEMPDIR",  ##not supported
        "-i ${files_in[$TASK_ID]}",
        "-o ${files_out[$TASK_ID]}"
      )
    ),
    arraysize = num_shards
  )  
}








################################################################################
################ Bascet command line tools: isolate genomes ####################
################################################################################


###############################################
####### Add raw FASTQs of isolates, enabling them to be treated as cells, assembled etc
# todo think about what to name this file. maybe separate from cell fastq to enable easy rerunning
# todo allow this function to be called multiple times?
BascetAddIsolateRawFastq <- function(bascetRoot, listFastqR1, listFastqR2, names, runner, bascet_instance=bascet_instance.default){
  
  
}




###############################################
####### Add isolate genomes, treating them as assembled, enabling clustering, comparison with cells, etc
# todo think about what to name this file. maybe separate from cell assembly to enable easy rerunning
# todo allow this function to be called multiple times?
BascetAddAssembledIsolate <- function(bascetRoot, listFasta, names, runner, bascet_instance=bascet_instance.default){
  
  
}




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
    showProgress=TRUE
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
    output[[cellname]] <- aggrFunction(bascetFile, cellname)
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
################ yet to classify ###############################################
################################################################################









