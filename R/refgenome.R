 

################################################################################
################## Alignment ###################################################
################################################################################


###############################################
#' Figure out if a BAM-file is a paired alignment or not
#' 
#' @param fname Name of BAM-file
#' @param bascetInstance Deprecated; ignored
#' @param maxRecords Maximum number of BAM records to scan
#' 
#' @return TRUE if the BAM-file is a paired alignment
#' @export
isBamPairedAlignment <- function(
    fname,
    bascetInstance=GetDefaultBascetInstance(),
    maxRecords=10000
){
  #Check arguments
  stopifnot(file.exists(fname))
  stopifnot(is.numeric(maxRecords))
  stopifnot(length(maxRecords) == 1)
  stopifnot(maxRecords > 0)
  
  bam <- Rsamtools::BamFile(fname, yieldSize = as.integer(maxRecords))
  open(bam)
  on.exit(close(bam), add = TRUE)
  bam_records <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = "flag")
  )
  flags <- bam_records[[1]]$flag
  any(bitwAnd(flags, 0x1) != 0)
}




###############################################
#' Index a genome using BWA-MEM2 such that it can be used for alignment
#'
#' TODO: could check if genome is indexed already
#'
#' @param genomeFile Name of FASTA file holding genome sequence
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#'
#' @export
BascetIndexGenomeBWAMEM2 <- function(
    genomeFile,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  #Check arguments
  stopifnot(is.existing.fasta(genomeFile))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))

  #if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    #Produce the script and run the job
    RunJob(
      runner = runner,
      jobname = "Zbwa_index",
      bascetInstance = bascetInstance,
      cmd = JobScript(
        steps = list(
          JobBascetCommand(bascetInstance, list(
            "exttool", "bwa-mem2", "index",
            genomeFile
          ))
        )
      ),
      arraysize = 1
    )
}





###############################################
#' Index a genome using Bowtie2 such that it can be used for alignment
#' 
#' @param genomeFile Name of FASTA file holding genome sequence
#' @param numThreads Number of threads to use. Default is the maximum, taken from runner settings
#' @param overwrite Force overwriting of existing files. The default is to do nothing if files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#'
#' @export
BascetIndexGenomeBowtie2 <- function(
    genomeFile,
    numThreads=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  #Set number of threads if not given
  if(is.null(numThreads)) {
    numThreads <- as.integer(runner@ncpu)
  }

  #Check arguments
  stopifnot(is.existing.fasta(genomeFile))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  
  indexName <- genomeFile #paste0(genomeFile,".btindex")
  
  if(bascetCheckOverwriteOutput(indexName, overwrite)) {
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Zbt2_build",
      bascetInstance = bascetInstance,
      cmd = JobScript(
        steps = list(
          JobCommand("bowtie2-build", list(
            JobArg("--threads", numThreads, sep = " "),
            genomeFile,
            indexName
          ), prepend = bascetInstance@prependCmd)
        )
      ),
      arraysize = 1
    )
  } else {
    new_no_job()
  }
}




###############################################
#' Index a genome using STAR such that it can be used for alignment.
#' Note that Bascet STAR (star-rs) automatically detects if GTF and FASTA are gzipped
#' 
#' @param fastaFile Name of FASTA file holding genome sequence
#' @param gtfFile GFF file holding genome annotation
#' @param outDir A directory in which to store the index. This directory will be created
#' @param numThreads The number of threads to use for the STAR index, for each runner. Default is the maximum, taken from runner settings
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetIndexGenomeSTAR <- function(  
    fastaFile, 
    gtfFile,
    outDir,
    numThreads=NULL,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  #Set number of threads if not given
  if(is.null(numThreads)) {
    numThreads <- as.integer(runner@ncpu)
  }
  
  #Check arguments
  stopifnot(is.existing.fasta(fastaFile))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  if(!file.exists(gtfFile)){
    stop("Could not find genome GTF file")
  }
  fastaFile <- normalizePath(fastaFile, winslash = "/", mustWork = TRUE)
  gtfFile <- normalizePath(gtfFile, winslash = "/", mustWork = TRUE)
  outDir <- normalizeOutputPath(outDir)
  
  #Exit early if already done
  if(file.exists(outDir)){
    print("STAR output directory already exists; assuming that genome is already built")
    return(new_no_job())
  }

    
  #Produce the script and run the job.
  #Note: overwrite not supported. It is too dangerous to implement
  RunJob(
    runner = runner,
    jobname = "Zstar_index",
    bascetInstance = bascetInstance,
    cmd = JobScript(
      steps = list(
        JobEnsureDir(outDir),
        JobEcho("Indexing genome with STAR..."),
        JobBascetCommand(bascetInstance, list(
          "exttool", "STAR",
          JobArg("--runThreadN", as.character(numThreads), sep = " "),
          JobArg("--runMode", "genomeGenerate", sep = " "),
          JobArg("--genomeDir", outDir, sep = " "),
          JobArg("--genomeFastaFiles", fastaFile, sep = " "),
          JobArg("--sjdbGTFfile", gtfFile, sep = " ")
        ))
      )
    ),
    arraysize = 1
  )
}




###############################################
#' Index a genome using minimap2 such that it can be used for alignment
#'
#' @param fastaFile Name of FASTA file holding genome sequence
#' @param indexFile Path to write the resulting .mmi index
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#'
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetIndexGenomeMinimap2 <- function(
    fastaFile,
    indexFile,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  #Check arguments
  stopifnot(is.existing.fasta(fastaFile))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))

  #Exit early if already done
  if(file.exists(indexFile)){
    print("minimap2 index already exists; skipping")
    return(new_no_job())
  }

  RunJob(
    runner = runner,
    jobname = "Zmm2_index",
    bascetInstance = bascetInstance,
    cmd = JobScript(
      steps = list(
        JobBascetCommand(bascetInstance, list(
          "exttool", "minimap2",
          JobArg("-d", indexFile, sep = " "),
          fastaFile
        ))
      )
    ),
    arraysize = 1
  )
}






###############################################
#' Check if a file is a FASTQ file
#' 
#' @param fname Path to file
#' 
#' @return TRUE if the file is some type of FASTQ
#' @noRd
isFastq <- function(
    fname
) {
  stringr::str_ends(fname, stringr::fixed("fq.gz")) ||  
    stringr::str_ends(fname, stringr::fixed("fastq.gz"))
}

###############################################
#' Check if a file is a paired FASTQ file.
#' Panics if the file is not a FASTQ at all
#' 
#' @param fname Path to file
#' 
#' @return TRUE if the file is a paired FASTQ
#' @noRd
isPairedFastq <- function(
    fname
) {
  if(isFastq(fname)){
    stringr::str_ends(fname, stringr::fixed("R1.fq.gz")) ||  
      stringr::str_ends(fname, stringr::fixed("R1.fastq.gz")) ||
      stringr::str_ends(fname, stringr::fixed("R2.fq.gz")) ||  
      stringr::str_ends(fname, stringr::fixed("R2.fastq.gz"))
  } else {
    stop("The file is not any type of FASTQ")
  }
}



###############################################
#' Get corresponding R2 file. Assumes that the input file is R1
#' 
#' @param fname Path to file
#' 
#' @return Path to R2 file
#' @noRd
getFastqR2fromR1 <- function(
    fname
) {

  dir <- dirname(fname)
  bname <- basename(fname)

  bname <- stringr::str_replace(bname, stringr::fixed("R1.fq.gz"), "R2.fq.gz")  
  bname <- stringr::str_replace(bname, stringr::fixed("R1.fastq.gz"), "R2.fastq.gz")  

  file.path(dir, bname)
}





###############################################
#' Generate a bigwig out of all reads in a sorted BAM. Note that the caller
#' is responsible for sorting the BAM first
#' 
#' This function is mainly for QC purposes.
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard (BAM-file)
#' @param outputName Name of output shard (BIGWIG-file)
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param numThreads Number of threads for BAM decompression and BigWig writing
#' @param totalMem Total memory to allocate
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A runner job (details depends on runner)
#' @export
BascetAlignmentToBigwig <- function(  
    bascetRoot, 
    inputName="aligned_pos", 
    outputName="pileup",
    overwrite=FALSE,
    numThreads=NULL,
    totalMem=NULL,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  if(is.null(numThreads)){
    numThreads <- as.integer(runner@ncpu)
  }

  #Check arguments 
  stopifnot(dir.exists(bascetRoot))
  bascetRoot <- normalizeBascetRoot(bascetRoot)
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.logical(overwrite))
  stopifnot(is.valid.threadcount(numThreads))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))

  #Check memory sizes
  totalMem <- checkTotalMemArg(totalMem, runner, bascetInstance)
  
  #Figure out input and output file names
  input_shards <- detectShardsForFile(bascetRoot, inputName) 
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards) 
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "bigwig", num_shards)
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    
    ### Build command
    cmd <- JobScript(
      vars = list(
        files_in = inputFiles,
        files_out = outputFiles
      ),
      steps = list(
        JobBascetCommand(bascetInstance, list(
          "tobigwig",
          JobArg("--in", JobVar("files_in"), sep = " "),
          JobArg("--out", JobVar("files_out"), sep = " "),
          JobMaybeArg("--memory", totalMem, format_size_bascet),
          JobArg("--threads", numThreads, sep = " ")
        ))
      )
    )

    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Ztobw",
      bascetInstance = bascetInstance,
      cmd = cmd,
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
}










###############################################
#' Filter an alignment (BAM-file).
#' 
#' This is typically used to either remove host DNA, or keep reads mapping to a known reference.
#' If the BAM-file has paired reads then BOTH reads need to be mapped (flag 0x2); otherwise (flag 0x4)
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param numThreads Number of threads to use. Default is the maximum, taken from runner settings
#' @param inputName Name of input shards (BAM-file format)
#' @param outputName Name of output shards (BAM-file format)
#' @param keepMapped Keep the mapped reads (TRUE) or unmapped (FALSE)
#' @param totalMem Total memory to allocate
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetFilterAlignment <- function(  
    bascetRoot, 
    numThreads=NULL,
    inputName, 
    outputName,
    keepMapped=FALSE,
    totalMem=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  #Set number of threads if not given
  if(is.null(numThreads)) {
    numThreads <- as.integer(runner@ncpu)
  }
  
  #Check input arguments 
  stopifnot(dir.exists(bascetRoot))
  bascetRoot <- normalizeBascetRoot(bascetRoot)
  stopifnot(is.valid.threadcount(numThreads))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.logical(keepMapped))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))

  #Check memory sizes
  totalMem <- checkTotalMemArg(totalMem, runner, bascetInstance)
  
  #Figure out input and output file names
  input_shards <- detectShardsForFile(bascetRoot, inputName) 
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards) 
  
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "bam", num_shards)

  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    
    keep_arg <- if(keepMapped) "mapped" else "unmapped"
    
    ### Build command
    cmd <- JobScript(
      vars = list(
        files_in = inputFiles,
        files_out = outputFiles
      ),
      steps = list(
        JobBascetCommand(bascetInstance, list(
          "filterbam",
          JobArg("--in", JobVar("files_in"), sep = " "),
          JobArg("--out", JobVar("files_out"), sep = " "),
          JobArg("--temp", JobEnv("BASCET_TEMPDIR"), sep = " "),
          JobArg("--threads", numThreads, sep = " "),
          JobMaybeArg("--memory", totalMem, format_size_bascet),
          JobArg("--keep", keep_arg, sep = " ")
        ))
      )
    )
      
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "ZfilterAln",
      bascetInstance = bascetInstance,
      cmd = cmd,
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
}







###############################################
#' Align from FASTQ, generate sorted and indexed BAM file
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param useReference Name of the BWA reference to use
#' @param numThreads Number of threads to use for each runner. Default is the maximum, taken from runner settings
#' @param totalMem Total memory to allocate
#' @param inputName Name of input shard
#' @param outputName Output shard name prefix. Defaults to `<outputName>_cell`
#'   for the unsorted/cell BAM and `<outputName>_pos` for the position-sorted BAM.
#' @param outputNameBAMcell Name of cell-sorted BAMs. If NULL, derived from `outputName`.
#' @param outputNameBAMpos Name of pos-sorted BAMs. If NULL, derived from `outputName`.
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param aligner Which aligner to use: "BWAMEM2", "STAR", or "minimap2"
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#'
#' @export
BascetAlignToReference <- function(
    bascetRoot, 
    useReference,
    numThreads=NULL,
    totalMem=NULL,
    inputName="filtered",
    outputName="aligned",
    outputNameBAMcell=NULL,
    outputNameBAMpos=NULL,
    overwrite=FALSE,
    aligner=c(NULL, "BWAMEM2", "STAR", "minimap2"),
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  #Set number of threads if not given
  if(is.null(numThreads)) {
    numThreads <- as.integer(runner@ncpu)
  }
  
  #Check input arguments 
  stopifnot(dir.exists(bascetRoot))
  bascetRoot <- normalizeBascetRoot(bascetRoot)
  stopifnot(is.valid.threadcount(numThreads))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  if(is.null(outputNameBAMcell)) {
    outputNameBAMcell <- paste0(outputName, "_cell")
  }
  if(is.null(outputNameBAMpos)) {
    outputNameBAMpos <- paste0(outputName, "_pos")
  }
  stopifnot(is.valid.shardname(outputNameBAMcell))
  stopifnot(is.valid.shardname(outputNameBAMpos))
  stopifnot(is.logical(overwrite))
  
  aligner <- match.arg(aligner)
  stopifnot("`aligner` must be specified." = !is.null(aligner))
  
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  
  
  #Figure out input and output file names  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards)
  outputFilesBAMunsorted  <- makeOutputShardNames(bascetRoot, outputNameBAMcell,  "bam", num_shards)
  outputFilesBAMsorted    <- makeOutputShardNames(bascetRoot, outputNameBAMpos, "bam", num_shards)
  
  if(aligner=="BWAMEM2"){
    if(!file.exists(useReference)){
      stop("BWAMEM2 reference file does not exist")
    }
  } else if (aligner=="minimap2") {
    if(!file.exists(useReference)){
      stop("minimap2 reference file does not exist")
    }
  } else if(aligner=="STAR") {
    #stop("not yet implemented")
    if(!dir.exists(useReference)) {  
      ### TODO could also check contents of this dir
      stop("STAR reference file does not exist")
    }
  } else {
    stop(paste("Unknown aligner: ", aligner))
  }
  
  #Check memory sizes
  totalMem <- checkTotalMemArg(totalMem, runner, bascetInstance)
  
  if(bascetCheckOverwriteOutput(outputFilesBAMsorted, overwrite)) {
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Zaln",
      bascetInstance = bascetInstance,
      cmd = JobScript(
        vars = list(
          files_in = inputFiles,
          files_out_unsorted = outputFilesBAMunsorted,
          files_out_sorted = outputFilesBAMsorted
        ),
        steps = list(
          if(!overwrite) JobSkipIfFileExists(JobVar("files_out_sorted")),
          JobBascetCommand(bascetInstance, list(
            "align",
            JobArg("--in", JobVar("files_in")),
            JobArg("--unsorted", JobVar("files_out_unsorted")),
            JobArg("--sorted", JobVar("files_out_sorted")),
            JobArg("--temp", JobEnv("BASCET_TEMPDIR")),
            JobArg("--genome", useReference),
            JobMaybeArg("--memory", totalMem, format_size_bascet),
            JobMaybeArg("--threads", numThreads),
            JobArg("--aligner", aligner)
          ))
        )
      ),
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
}






###############################################
#' Take aligned BAM file and produce Fragments.tsv.gz, compatible with Signac ATAC-seq style analysis
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A job, producing a type of Fragments.tsv.gz
#' @export
BascetBam2Fragments <- function(
    bascetRoot, 
    inputName="aligned_pos",
    outputName="fragments", 
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  #Check input arguments 
  stopifnot(dir.exists(bascetRoot))
  bascetRoot <- normalizeBascetRoot(bascetRoot)
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  
  #Get input names  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  if(length(input_shards)==0){
    stop("Found no input files")
  }
  num_shards <- length(input_shards)
  inputFiles <- file.path(bascetRoot, input_shards)
  
  #One output per input pair of reads  
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "tsv.gz", num_shards)  
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    RunJob(
      runner = runner, 
      jobname = "Zbam2frag",
      bascetInstance = bascetInstance,
      cmd = JobScript(
        vars = list(
          files_in = inputFiles,
          files_out = outputFiles
        ),
        steps = list(
          if(!overwrite) JobSkipIfFileExists(JobVar("files_out")),
          JobBascetCommand(bascetInstance, list(
            "bam2fragments",
            JobArg("-t", JobEnv("BASCET_TEMPDIR")),
            JobArg("-i", JobVar("files_in")),
            JobArg("-o", JobVar("files_out"))
          ))
        )
      ),
      arraysize = length(inputFiles)
    )
  } else {
    new_no_job()
  }
}





###############################################
#' From aligned BAM file, compute counts per chromosome
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param minMatching Disregard reads having fewer than specified matches, based on CIGAR string
#' @param removeDuplicates Deduplicate reads
#' @param removeMultimapper Remove reads for a cell if they map to multiple places
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetCountChrom <- function(
    bascetRoot,
    inputName="aligned_pos",
    outputName="chromcount", 
    minMatching=0,
    removeDuplicates=TRUE,
    removeMultimapper=TRUE,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  #Check input arguments 
  stopifnot(dir.exists(bascetRoot))
  bascetRoot <- normalizeBascetRoot(bascetRoot)
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.numeric(minMatching))
  stopifnot(is.logical(removeDuplicates))
  stopifnot(is.logical(removeMultimapper))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  
  #Get input names  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  if(length(input_shards)==0){
    stop("Found no input files")
  }
  num_shards <- length(input_shards)
  inputFiles <- file.path(bascetRoot, input_shards)
  
  #One output per input alignment
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "h5", num_shards)  
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    RunJob(
      runner = runner, 
      jobname = "Zcountchrom",
      bascetInstance = bascetInstance,
      cmd = JobScript(
        vars = list(
          files_in = inputFiles,
          files_out = outputFiles
        ),
        steps = list(
          if(!overwrite) JobSkipIfFileExists(JobVar("files_out")),
          JobBascetCommand(bascetInstance, list(
            "countchrom",
            JobArg("--min-matching", minMatching),
            if(removeDuplicates) JobArg("--remove-duplicates"),
            if(removeMultimapper) JobArg("--remove-multimapper"),
            JobArg("-t", JobEnv("BASCET_TEMPDIR")),
            JobArg("-i", JobVar("files_in")),
            JobArg("-o", JobVar("files_out"))
          ))
        )
      ),
      arraysize = length(inputFiles)
    )
  } else {
    new_no_job()
  }
}





###############################################
#' Using Tabix, get list of sequences in a fragment file
#' 
#' @param fragpath List to fragment file
#' @param bascetInstance A Bascet instance
#' 
#' @return The raw result of tabix
#' @export
TabixGetFragmentsSeqs <- function(
    fragpath,
    bascetInstance=GetDefaultBascetInstance()
) {
  #Check input arguments 
  stopifnot(file.exists(fragpath))
  stopifnot(is.bascet.instance(bascetInstance))
  
  system(paste(bascetInstance@prependCmd,"tabix -l ", fragpath), intern = TRUE)
}


###############################################
#' From a fragments file, get a chromatin assay for Signac.
#' 
#' Note: signac is dirt slow at counting as of writing. It might be scalable enough
#' for certain inputs and tasks, but we still provide the option of using it.
#' 
#' @param fragpath List to fragment file
#' 
#' @return A ChromatinAssay
#' @export
FragmentsToSignac <- function(
    fragpath
) {
  #Check input arguments 
  stopifnot(file.exists(fragpath))

  #### Index if needed; this should already have been done but added here just in case
  fragpath_index <- paste(fragpath,".tbi",sep="")
  if(!file.exists(fragpath_index)){
    print("Indexing fragment file")
    system(paste("tabix -p bed ",fragpath))  ########### TODO: use singularity container?
  }
  
  #### Count fragments / see which cells are present. Possible to add a cutoff here if we wish
  counts <- CountFragments(fragments = fragpath)  
  colnames(counts)[1] <- "cb"
  counts <- counts[order(counts$reads_count, decreasing = TRUE),]
  
  #knee plot analysis
  if(FALSE){
    kneedf <- data.frame(cnt=sort(counts$reads_count, decreasing = TRUE))
    kneedf$index <- 1:nrow(kneedf)
    ggplot(kneedf, aes(index, cnt)) + 
      geom_line() + 
      scale_x_log10() + 
      scale_y_log10()
  }
  
  keep_cells <- counts$cb
  
  #### Create a dummy assay
  stupidmat <- matrix(0, nrow = 2, ncol=length(keep_cells))
  rownames(stupidmat) <- c("chr1:1-100", "chr1:200-300")
  colnames(stupidmat) <- keep_cells
  
  chromAssay <- CreateChromatinAssay(
    counts = as.sparse(stupidmat),
    sep = c(":", "-"),
    fragments = fragpath,
    min.cells = 0,
    min.features = 0
  )  
  chromAssay
}



###############################################
#' From a Signac chromatin assay with fragments, for each cell, count how many reads per chromosome
#' 
#' @param chromAssay A ChromatinAssay
#' @param bascetInstance A Bascet instance
#' 
#' @return a FeatureMatrix 
#' @export
FragmentCountsPerChrom <- function(
    chromAssay,
    bascetInstance=GetDefaultBascetInstance()
){
  
  #Figure out where the fragment file is
  fr <- Fragments(chromAssay)[[1]]
  
  #Get the name of chromosomes
  all_seqid <- TabixGetFragmentsSeqs(fr@path, bascetInstance)
  grange <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqid=all_seqid,
    name=all_seqid,
    start=1,
    end=536870912  #for grange this must be <2**31. for scanTabix, this must be <=536870912
  ))
  
  # Quantify counts over each chromosome
  reg_counts <- FeatureMatrix(
    fragments = Fragments(chromAssay),
    features = grange,
    cells = colnames(fr@cells)
  )  
  
  reg_counts
}




###############################################
#' From a Signac chromatin assay with fragments, for each cell, count how many reads per chromosome.
#' This function directly returns an assay that can be added to a Seurat multimodal object
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' 
#' @return A seurat object holding the counts
#' @export
FragmentCountsPerChromAssay <- function(
    bascetRoot,
    inputName="fragments.1.tsv.gz" ### TODO, detect, merge all of them
) {
  #Check input arguments 
  stopifnot(dir.exists(bascetRoot))
  bascetRoot <- normalizeBascetRoot(bascetRoot)

  stopifnot(is.character(inputName), length(inputName) == 1)  ##### need a better test in the future!
  stopifnot(!is.na(inputName), nzchar(inputName)) ##### need a better test in the future!
  
  fragpath <- file.path(bascetRoot,inputName)
  stopifnot(file.exists(fragpath))

  chromAssay <- FragmentsToSignac(fragpath)
  chromAssay <- FragmentCountsPerChrom(chromAssay)
  
  chromAssay_hack <- CreateAssay5Object(counts = chromAssay)
  rownames(chromAssay_hack) <- stringr::str_split_i(rownames(chromAssay_hack),"-1-",1)  #evil hack, fragile

  chromAssay_hack
}




###############################################
#' Produce a count matrix on strain level
#' 
#' This function will likely be depreciated in the future, as it is a bit too specific to keep in the library
#' 
#' @param adata A Seurat object
#' @param mapSeq2strain A mapping from sequence to strain
#' 
#' @return ...
#' @export
ChromToSpeciesCount <- function(
    adata, 
    mapSeq2strain
){
  #Check input arguments 
  stopifnot(inherits(adata, "Seurat"))
  stopifnot(is.data.frame(mapSeq2strain))
  stopifnot(all(c("id", "strain") %in% colnames(mapSeq2strain)))
  
  mat_cnt <- adata@assays[[DefaultAssay(adata)]]$counts
  stopifnot(!is.null(mat_cnt))
  stopifnot(nrow(mat_cnt) > 0, ncol(mat_cnt) > 0)
  stopifnot(any(mapSeq2strain$id %in% rownames(mat_cnt)))
  
  unique_strains <- unique(mapSeq2strain$strain)
  strain_cnt <- matrix(NA, ncol=ncol(mat_cnt), nrow=length(unique_strains))
  rownames(strain_cnt) <- unique_strains
  colnames(strain_cnt) <- colnames(mat_cnt)
  for(i in 1:nrow(strain_cnt)){
    cur_cols <- mapSeq2strain$id[mapSeq2strain$strain %in% rownames(strain_cnt)[i]]
    print(cur_cols)
    strain_cnt[i,] <- colSums(mat_cnt[rownames(mat_cnt) %in% cur_cols,,drop=FALSE])
  }
  
  CreateAssay5Object(counts = strain_cnt)
}






################################################################################
############ RNA-seq style feature counting from fragments.tsv #################
################################################################################


###############################################
#' Obtain a feature matrix (as seurat object) given an seurat object having Fragments associated
#' 
#' @param adata Seurat object, with FeatureMatrix
#' @param grangeGene A grange object telling where to count
#' @return A Seurat AssayObject
#' 
#' @export
CountGrangeFeatures <- function(
    adata, 
    grangeGene
){
  #Check input arguments 
  stopifnot(inherits(adata, "Seurat"))
  stopifnot(length(Fragments(adata)) > 0)
  stopifnot(inherits(grangeGene, "GRanges"))
  stopifnot("Name" %in% colnames(S4Vectors::mcols(grangeGene)))
  
  gene_counts <- FeatureMatrix(
    fragments = Fragments(adata),
    features = grangeGene,
    cells = colnames(adata)
  )
  rownames(gene_counts) <- grangeGene$Name ## paste(grangeGene$gene_biotype,grangeGene$Name) #Name in an easy manner
  
  # create a new assay using the MACS2 peak set 
  CreateAssayObject(
    counts = gene_counts[!duplicated(rownames(gene_counts)),],
    min.cells = 0,  #filter later!!!
    min.features = 0
  )
}













###############################################
#' From aligned BAM file, compute counts per feature
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param gffFile GFF-like file to use for feature annotation
#' @param useFeature Feature type in the file to count
#' @param attrGeneId Attribute field to use this on for gene ID
#' @param attrGeneName Attribute field to use this on for gene name
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A job, executing the counting 
#' @export
BascetCountFeature <- function(
    bascetRoot,
    inputName="aligned_pos",
    outputName="featurecount", 
    gffFile,
    useFeature="gene",
    attrGeneId="gene_id",  #"ID" for yersinia
    attrGeneName="name",
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  #Check input arguments 
  stopifnot(dir.exists(bascetRoot))
  bascetRoot <- normalizeBascetRoot(bascetRoot)
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(file.exists(gffFile))
  stopifnot(is.character(useFeature))
  stopifnot(is.character(attrGeneId))
  stopifnot(is.character(attrGeneName))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  
  #Get input name
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  if(length(input_shards)==0){
    stop("Found no input files")
  }
  num_shards <- length(input_shards)
  inputFiles <- file.path(bascetRoot, input_shards)
  
  #One output per input alignment
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "h5", num_shards)  
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    RunJob(
      runner = runner, 
      jobname = "Zcountf",
      bascetInstance = bascetInstance,
      cmd = JobScript(
        vars = list(
          files_in = inputFiles,
          files_out = outputFiles
        ),
        steps = list(
          if(!overwrite) JobSkipIfFileExists(JobVar("files_out")),
          JobBascetCommand(bascetInstance, list(
            "countfeature",
            JobArg("-g", gffFile),
            JobArg("--use-feature", useFeature),
            JobArg("--attr-id", attrGeneId),
            JobArg("--attr-name", attrGeneName),
            JobArg("-t", JobEnv("BASCET_TEMPDIR")),
            JobArg("-i", JobVar("files_in")),
            JobArg("-o", JobVar("files_out"))
          ))
        )
      ),
      arraysize = length(inputFiles)
    )
  } else {
    new_no_job()
  }
}







################################################################################
################## SNP analysis ################################################
################################################################################

#https://academic.oup.com/bioinformatics/article/37/23/4569/6272512

#like KRAKEN, wrap cellsnp-lite

# is b needed?
# cellsnp-lite -s $BAM -b $BARCODE -O $OUT_DIR -p 10 --minMAF 0.1 --minCOUNT 100 --gzip


# singularity run ~/github/bascet/singularity/bascet.sif cellsnp-lite -s mutans_aligned.1.bam -O cellsnp.1.out -p 1 --minMAF 0.1 --minCOUNT 100 --gzip -b bclist.csv
#C1_A5_C8_A11
#F3_G6_H8_C11


# singularity run ~/mystore/bascet.sif cellsnp-lite -s mutans_aligned.1.bam -O cellsnp.1.out -p 1 --minMAF 0.1 --minCOUNT 100 --gzip -b bclist.csv
# creates output dir



###############################################
#' Align from FASTQ, generate sorted and indexed BAM file
#' 
#' TODO Should be managed by Bascet mapshard system, with automatic input conversion. unaligned file should be made temp and removed
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param numThreads Number of threads to use for each runner. Default is the maximum, taken from runner settings
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @export
BascetRunCellSNP <- function(  
    bascetRoot, 
    inputName="aligned_pos", 
    outputName="cellsnp", 
    numThreads=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  #Set number of threads if not given
  if(is.null(numThreads)) {
    numThreads <- as.integer(runner@ncpu)
  }
  
  #Check input arguments 
  stopifnot(dir.exists(bascetRoot))
  bascetRoot <- normalizeBascetRoot(bascetRoot)
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
  
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "out", num_shards)
  
  listcellFiles <- makeOutputShardNames(bascetRoot, outputName, "listcell", num_shards)
  listchromFiles <- makeOutputShardNames(bascetRoot, outputName, "listchrom", num_shards)

  ### Build command: basic alignment
  cmd <- JobScript(
    vars = list(
      files_in = inputFiles,
      files_out = outputFiles,
      files_listcell = listcellFiles,
      files_listchrom = listchromFiles
    ),
    steps = list(
      if(!overwrite) JobSkipIfFileExists(JobVar("files_listcell")),
      JobRedirect(
        JobPipeline(
          JobCommand("samtools", list("idxstats", JobVar("files_in")), prepend = bascetInstance@prependCmd),
          JobCommand("head", list("-n", "-1")),
          JobCommand("cut", list("-f", "1"))
        ),
        stdout = JobVar("files_listchrom")
      ),
      JobRedirect(
        JobPipeline(
          JobCommand("samtools", list("view", JobVar("files_in")), prepend = bascetInstance@prependCmd),
          JobCommand("sed", list("-e", "s/^.*CB:Z://")),
          JobCommand("sed", list("-e", "s/\t.*//")),
          JobCommand("sort"),
          JobCommand("uniq")
        ),
        stdout = JobVar("files_listcell")
      ),
      JobCommand("cellsnp-lite", list(
        JobArg("-s", JobVar("files_in"), sep = " "),
        JobArg("-p", numThreads, sep = " "),
        JobArg("--genotype"),
        JobArg("--chrom", JobFileLinesCsv(JobVar("files_listchrom")), sep = " "),
        JobArg("--gzip"),
        JobArg("-b", JobVar("files_listcell"), sep = " "),
        JobArg("-O", JobVar("files_out"), sep = " ")
      ), prepend = bascetInstance@prependCmd)
    )
  )
  
  
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Zcellsnp",
      bascetInstance = bascetInstance,
      cmd = cmd,
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
}






###############################################
### Helper function
ReadCellSNPmatrix_one <- function(
    basedir
) {
  #Check input arguments 
  stopifnot(dir.exists(basedir))
  
  #cellSNP.base.vcf  one line per feature
  snp_info <- read.table(file.path(basedir,"cellSNP.base.vcf.gz"))
  colnames(snp_info) <- c("chrom","pos","id","ref","alt","qual","filter","info")
  rownames(snp_info) <- sprintf(
    "%s_%s_%s_to_%s",
    snp_info$chrom, snp_info$pos, snp_info$ref, snp_info$alt
  )
  
  snp_cellid <- readLines(file.path(basedir,"cellSNP.samples.tsv"))
  length(snp_cellid)
  
  #depth of ALT alleles
  snp_ad <- Matrix::readMM(file.path(basedir,"cellSNP.tag.AD.mtx"))
  colnames(snp_ad) <- snp_cellid
  rownames(snp_ad) <- rownames(snp_info)
  
  #depth of ALT+REF alleles
  snp_dp <- Matrix::readMM(file.path(basedir,"cellSNP.tag.DP.mtx"))
  colnames(snp_dp) <- snp_cellid
  rownames(snp_dp) <- rownames(snp_info)
  
  #set feature names
  list(
    snp_ad=snp_ad,
    snp_dp=snp_dp,
    snp_info=snp_info
  )
}


###############################################
#' Read a count matrix as produced by CellSNP, but as shards
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param listCells Optional list of cells to extract
#' @param verbose Show process status for debugging purposes
#' 
#' @return Count matrix as sparseMatrix
#' @export
ReadCellSNPmatrix <- function(
    bascetRoot, 
    inputName,
    listCells=NULL,
    verbose=FALSE
){
  #Check input arguments 
  stopifnot(dir.exists(bascetRoot))
  bascetRoot <- normalizeBascetRoot(bascetRoot)
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.listcells(listCells))
  stopifnot(is.logical(verbose))

  print("Loading CellSNP file")
  
  #bascetRoot <- "/husky/henriksson/atrandi/v4_wgs_saliva1"
  #inputName <- "cellsnp"
  
  #Figure out input file names  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards)
  #if(tools::file_ext(inputFiles[1])!="h5"){
  #  stop("Wrong input format. should be hd5")
  #}
  
  #Load individual matrices. Sizes may not match
  if(verbose){
    print("Loading matrices")
  }
  list_mat <- list()
  for(f in inputFiles){
    list_one <- ReadCellSNPmatrix_one(f) 
    
    mat <- Matrix::t(list_one$snp_ad)
    
    if(!is.null(listCells)){
      rownames(mat)
      mat <- mat[rownames(mat) %in% listCells,,drop=FALSE]
    }
    
    if(verbose){
      print(dim(mat))
    }
    list_mat[[f]] <- mat
  }
  
  #Find union of features. Assume that no filtering was done such that we can assume 0
  #count for a feature not present in another shard
  if(verbose){
    print("Union of features across matrices")
  }
  all_colnames <- sort(unique(unlist(lapply(list_mat, colnames))))
  if(verbose){
    print(all_colnames)
  }
  num_col <- length(all_colnames)
  map_name_to_i <- data.frame(row.names = all_colnames, ind=1:length(all_colnames))
  if(verbose){
    print(map_name_to_i)
  }
  
  #Make sizes compatible
  if(verbose){
    print("Concatenating")
  }
  list_resized_mat <- list()
  for(f in inputFiles){
    if(verbose){
      print(f)
    }
    mat <- list_mat[[f]]
    new_mat <- MatrixExtra::emptySparse(nrow = nrow(mat), ncol = num_col, format = "R", dtype = "d")
    new_mat[1:nrow(mat), map_name_to_i[colnames(mat),]] <- MatrixExtra::as.csr.matrix(mat)  #manually look up column names!  #here, x[.,.] <- val : x being coerced from Tsparse* to CsparseMatrix
    rownames(new_mat) <- rownames(mat)
    colnames(new_mat) <- all_colnames
    # print(dim(new_mat))
    list_resized_mat[[f]] <- new_mat
  }
  
  #Concatenate matrices
  allmat <- do.call(rbind, list_resized_mat) #TODO check that above worked properly!
  
  allmat
}
