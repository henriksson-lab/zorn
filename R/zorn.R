




################################################################################
################ Bascet command line tool: debarcoding raw fastq ###############
################################################################################



###############################################
#' Detect metadata for raw input FASTQ files
#' 
#' @param rawRoot Path to folder with FASTQ files
#' @returns A data frame with metadata for the raw input files
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
BascetGetRawAtrandiWGS <- function(
    bascetRoot, 
    rawmeta, 
    outputName="debarcoded", 
    outputNameIncomplete="incomplete_reads", 
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
    withdata = c(
    ),
    cmd = c(
      shellscript_make_bash_array("files_r1",file.path(rawmeta$dir, rawmeta$r1)),
      shellscript_make_bash_array("files_r2",file.path(rawmeta$dir, rawmeta$r2)),
      shellscript_make_bash_array("files_out",outputFilesComplete),
      shellscript_make_bash_array("files_out_incomplete",outputFilesIncomplete),
      paste(
        bascet_instance@bin, 
        "getraw",
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
    withdata = c(), #drop later? 
    cmd = c(
      shellscript_make_files_expander("CELLFILE", list_cell_for_shard),
      shellscript_make_bash_array("files_out", outputFiles),
      paste(
        bascet_instance@bin,
        "shardify", 
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
#' Generate zip file with fastq from each cell
# this command should be able to take multiple debarcoded files as input
BascetPartition <- function(bascetRoot, inputName="debarcoded", outputName="rawreads", runner, bascet_instance=bascet_instance.default){
  
  #Figure out input and output file names  
  inputFiles <- detect_shards_for_file(bascetRoot, inputName)
  num_shards <- length(inputFiles)
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "zip", num_shards)
  
  #Run the job
  RunJob(
    runner = runner, 
    jobname = "bascet_partition",
    withdata = c(
      shellscript_make_bash_array("files_in",inputFiles),
      shellscript_make_bash_array("files_out",outputFiles) ####### drop!
    ),
    cmd = paste(bascet_instance@bin, "partitionTODO -i $files_in[$TASK_ID] -o $files_out[$TASK_ID]"),
    arraysize = num_shards
  )
  
}


###############################################
#' Run Spades to assemble the genomes
BascetAssemble <- function(bascetRoot, inputName="rawreads", outputName="assembled", runner, bascet_instance=bascet_instance.default){
  
  #Figure out input and output file names  
  inputFiles <- detect_shards_for_file(bascetRoot, inputName)
  num_shards <- length(inputFiles)
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "zip", num_shards)
  
  #Run the job
  RunJob(
    runner = runner, 
    jobname = "bascet_asm",
    withdata = c(
      shellscript_make_bash_array("files_in",inputFiles),
      shellscript_make_bash_array("files_out",outputFiles)
    ),
    cmd = paste("bascet assemble TODO --in files_in[$TASK_ID] --o files_out[$TASK_ID]"),
    arraysize = num_shards
  )
  
}


###############################################
#' Build kmer database
BascetCount <- function(bascetRoot, inputName="assembled", outputName="kmers", runner, bascet_instance=bascet_instance.default){
  
  
  #Figure out input and output file names  
  inputFiles <- detect_shards_for_file(bascetRoot, inputName)
  num_shards <- length(inputFiles)
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "zip", num_shards)
  
  #Run the job
  RunJob(
    runner = runner, 
    jobname = "bascet_kmc",
    withdata = c(
      shellscript_make_bash_array("files_in",inputFiles),
      shellscript_make_bash_array("files_out",outputFiles)
    ),
    cmd = paste("bascet count TODO --in files_in[$TASK_ID] --o files_out[$TASK_ID]"),
    arraysize = num_shards
  )
  
}



###############################################
#' Select kmers that appear useful for clustering
BascetFeaturise <- function(bascetRoot, inputName="kmers", outputName="features", runner, bascet_instance=bascet_instance.default){
  
  #Figure out input and output file names  
  inputFiles <- detect_shards_for_file(bascetRoot, inputName)
  num_shards <- length(inputFiles)
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "zip", num_shards)
  
  #Run the job
  RunJob(
    runner = runner, 
    jobname = "bascet_featurise",
    withdata = c(
      make_bash_array("files_in",inputFiles),
      make_bash_array("files_out",outputFiles)
    ),
    cmd = paste("bascet featurise TODO --in files_in[$TASK_ID] --o files_out[$TASK_ID]"),
    arraysize = num_shards
  )
  
}


###############################################
#' Build count table from kmer table and selected kmers
BascetQuery <- function(bascetRoot, inputKMER="kmers", inputFeatures="features", outputName="query", runner, bascet_instance=bascet_instance.default){
  
  #Figure out input and output file names  
  inputKMER <- detect_shards_for_file(bascetRoot, inputKMER)
  num_shards <- length(inputKMER)
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "zip", num_shards)
  
  #Note: Need to give all input files for each job
  inputFeatures_comma <- stringr::str_flatten(inputFeatures, collapse = ",")
  
  #Run the job
  RunJob(
    runner = runner, 
    jobname = "bascet_query",
    withdata = c(
      make_bash_array("files_kmer",inputKMER),
      make_bash_array("files_out",outputFiles)
    ),
    cmd = paste(bascet_instance@bin, "query --features ",inputFiles_comma," --kmers files_kmer[$TASK_ID] --o files_out[$TASK_ID]"), 
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
BascetMapCell <- function(
    bascetRoot, 
    withfunction, 
    inputName, 
    outputName, 
    runner, 
    bascet_instance=bascet_instance.default){
  
  #Figure out input and output file names  
  inputFiles <- detect_shards_for_file(bascetRoot, inputName)
  num_shards <- length(inputFiles)
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "zip", num_shards)
  
  #Run the job
  RunJob(
    runner = runner, 
    jobname = paste0("bascet_map_",withfunction),
    withdata = c(
      make_bash_array("files_in",inputFiles),
      make_bash_array("files_out",outputFiles)
    ),
    cmd = paste(bascet_instance@bin, "map -i $files_in[$TASK_ID] -o $files_out[$TASK_ID] -s ", withfunction),
    arraysize = num_shards
  )  
}




###############################################
#' Convenience function; alternative is to somehow implement as.data.frame
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
BascetAggregateMap <- function(
    bascetRoot, 
    bascetName, 
    aggrFunction){
  
  #Get file coordinates of all objects in zip file
  cellname_coord <- BascetCellNames(bascetRoot, bascetName)
  
  #Open the file, prep for reading
  bascetFile <- OpenBascet(bascetRoot, bascetName)
  
  #Loop over all files in the bascet
  output <- list()
  for(cellname in bascet_file@cellmeta$cell){
    output[[cellname]] <- aggrFunction(bascetFile, cellname)
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
################ yet to classify #####################################
################################################################################






library(ggplot2)

keep_cell_min_reads <- 100





if(FALSE){
  
  #tmp <- BascetReadFile(bascetFile, cellID, "out.csv", as="tempfile")  from aggr function
  
  
  one_bascet <- OpenBascet("/home/mahogny/jupyter/bascet/zorn/try_unzip","quast")
  thefile <- BascetReadFile(one_bascet, "a", "report.txt")
  thefile
  dat <- readLines(thefile)
  
  
  aggr.quast
  
  thefile <- BascetReadFile("/home/mahogny/jupyter/zorn/", "0-1-2-3", "out.txt")
  
  (BascetAggregateMap("/home/mahogny/jupyter/bascet/zorn/try_unzip","quast",aggr.quast))
  
  MapListAsDataFrame(BascetAggregateMap("/home/mahogny/jupyter/bascet/zorn/try_unzip","quast",aggr.quast))
  #both are from a!
  
}


#transposed_report.tsv




