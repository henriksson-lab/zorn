


###############################################
#' 
#' Run KMC3 for FASTQ files. KMC3 does not support named pipes, so FASTQ conversion is unavoidable
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param numThreads Number of threads to use for each runner. Default is the maximum, taken from runner settings
#' 
#' TODO
#' 
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @export
BascetRunKmcPerShard <- function(  
    bascetRoot, 
    inputName="asfq", 
    outputName="kmcpershard", 
    numThreads=NULL,
    overwrite=FALSE,
    kmerSize=31,
    totalMem=NULL,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  #Set number of threads if not given
  if(is.null(numThreads)) {
    numThreads <- as.integer(runner@ncpu)
  }
  
  #Check input arguments 
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.valid.threadcount(numThreads))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  
  
  #Check memory size
  if(!is.null(totalMem)) {
    totalMem <- parse_size_string(totalMem)
    stopifnot(totalMem > fs::fs_bytes("1Gb"))
  } else {
    #Take memory from runner if possible
    if(runner@mem!="") {
      totalMem <- parse_size_string(runner@mem) - fs::fs_bytes(bascetInstance@containerMem)
      stopifnot(totalMem > fs::fs_bytes("1Gb"))
    } else {
      print("Warning: Total memory was not specified. We strongly encourage doing this to ensure performance")
    }
  }
  totalMemGb <- as.integer(round(totalMem/1e6)) #KMC want memory in Gb
  stopifnot(totalMemGb>0)
  
  
  #Figure out input and output file names  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFilesR1 <- file.path(bascetRoot, input_shards) 
  inputFilesR2 <- getFastqR2fromR1(inputFilesR1)
  
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, NULL, num_shards)
  outputFilesPre <- makeOutputShardNames(bascetRoot, outputName, "kmc_pre", num_shards)
  
  ### Build command: basic alignment
  cmd <- c(
    shellscriptMakeBashArray("files_in_r1", inputFilesR1),
    shellscriptMakeBashArray("files_in_r2", inputFilesR2),
    shellscriptMakeBashArray("files_out", outputFiles),
    shellscriptMakeBashArray("files_out_pre", outputFilesPre),
    
    ### Abort early if needed    
    if(!overwrite) shellscriptCancelJobIfFileExists("${outputFilesPre[$TASK_ID]}"),

    ### Make list of input files for 
    paste(
      bascetInstance@prependCmd,
      "bash -c \"",
      "echo ${files_in_r1[$TASK_ID]} > ${BASCET_TEMPDIR}/listinput.txt; echo ${files_in_r2[$TASK_ID]} >> ${BASCET_TEMPDIR}/listinput.txt",
      "\""
    ),
    
    ### Run KMC
    paste(
      bascetInstance@prependCmd,
      "kmc",
      paste0("-k",kmerSize),
      paste0("-t",numThreads),
      paste0("-m",totalMemGb),
      "-ci6",
      "@${BASCET_TEMPDIR}/listinput.txt",
      "${files_out[$TASK_ID]}",
      "${BASCET_TEMPDIR}"
    )
  )

  #print(cmd)
  
  if(bascetCheckOverwriteOutput(outputFilesPre, overwrite)) { #### TODO this will not work; need to add suffix to files
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Zkmc",
      bascetInstance = bascetInstance,
      cmd = cmd,
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
}





###############################################
#' 
#' TODO!!!!!!!!!!!! need to handle MERGE shards
#' 
#' 
#' Merge KMC3 outputs into one
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param numThreads Number of threads to use for each runner. Default is the maximum, taken from runner settings
#' 
#' TODO
#' 
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @export
BascetMergeKmcPerShard <- function(  
    bascetRoot, 
    inputName="kmcpershard", 
    outputName="merge_kmcpershard", #no. should be same as input 
    numThreads=NULL,
    overwrite=FALSE,
    totalMem=NULL,
    kmcCi=6,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  #Set number of threads if not given
  if(is.null(numThreads)) {
    numThreads <- as.integer(runner@ncpu)
  }
  
  #Check input arguments 
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.valid.threadcount(numThreads))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  
  
  
  #Figure out input and output file names  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards)
  
  print(inputFiles)
  stop(666)
  inputFiles <- stringr::str_remove(inputFiles, stringr::fixed(".kmc_pre"),"")

  outputFiles <- makeOutputShardNames(bascetRoot, outputName, NULL, 1)
  outputFilesPre <- makeOutputShardNames(bascetRoot, outputName, "kmc_pre", 1) ### TODO: should be a merge shard
  
  ### Build KMC op-file
  index_inputs <- 1:length(inputFiles)
  opfile <- c(
    "INPUT:",
    sprintf("in%s = ", index_inputs, inputFiles),
    "OUTPUT:",
    paste("kmc_tot =", stringr::str_flatten(sprintf("in%s", index_inputs)," + "))
    )
  
  ######### TODO - need to create this file
  
  ### Build command
  cmd <- c(
    shellscriptMakeBashArray("files_in", inputFiles),
    shellscriptMakeBashArray("files_out", outputFiles),
    shellscriptMakeBashArray("files_out_pre", outputFilesPre),
    
    ### Abort early if needed    
    if(!overwrite) shellscriptCancelJobIfFileExists("${outputFilesPre[$TASK_ID]}"),
    
    ### TODO store op-file ########################################

    ### Run KMC
    paste(
      bascetInstance@prependCmd,
      "kmc_tools",
      paste0("-k",kmerSize),
      paste0("-t",numThreads),
      paste0("-m",totalMemGb),
      paste0("-ci",kmcCi),
      "@${BASCET_TEMPDIR}/listinput.txt",
      "${files_out[$TASK_ID]}",
      "${BASCET_TEMPDIR}"
    )
  )
  
  #print(cmd)
  
  if(bascetCheckOverwriteOutput(outputFilesPre, overwrite)) { #### TODO this will not work; need to add suffix to files
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Zkmc",
      bascetInstance = bascetInstance,
      cmd = cmd,
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
}

# 
# INPUT:
#   in1 = kmcpershard.1
# in2 = kmcpershard.2
# in3 = kmcpershard.3
# in4 = kmcpershard.4
# in5 = kmcpershard.5
# in6 = kmcpershard.6
# in7 = kmcpershard.7
# in8 = kmcpershard.8
# in9 = kmcpershard.9
# in10 = kmcpershard.10
# in11 = kmcpershard.11
# in12 = kmcpershard.12
# in13 = kmcpershard.13
# in14 = kmcpershard.14
# in15 = kmcpershard.15
# in16 = kmcpershard.16
# in17 = kmcpershard.17
# in18 = kmcpershard.18
# in19 = kmcpershard.19
# in20 = kmcpershard.20
# OUTPUT:
#   kmc_tot = in1 + in2 + in3 + in4 + in3 + in4 + in5 + in6 + in7 + in8 + in9 + in10 + in11 + in12 + in13 + in14 + in15 + in16 + in17 + in18 + in19 + in20
# 
# 
# 







