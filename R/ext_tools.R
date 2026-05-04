
################################################################################
################ SKESA #########################################################
################################################################################


###############################################
#' Run SKESA on reads of all cells.
#' This is a thin wrapper around BascetMapCell
#'
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param ... Additional arguments passed to \code{\link{BascetMapCell}}
#'
#' @return A job to be executed, or being executed, depending on runner settings
#' @seealso \code{\link{BascetMapCell}}
#' @export
BascetMapCellSKESA <- function(
    bascetRoot,
    inputName="filtered",
    outputName="contigs",
    ...
){
  BascetMapCell(
    bascetRoot=bascetRoot,
    withfunction="_skesa",
    inputName=inputName,
    outputName=outputName,
    ...
  )
}

###############################################
#' Run integrated SKESA on reads of all cells.
#'
#' This uses the bascet integrated skesa command rather than the old mapcell
#' system.
#'
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param numThreads Total thread budget. Defaults to the runner CPU count
#' @param numSkesaWorkers Number of cells to assemble concurrently
#' @param numSkesaCores Number of cores to give each SKESA assembly
#' @param numThreadsRead Threads used by the TIRP reader. If NULL, use the CLI default
#' @param totalMem Total memory to allocate
#' @param kmer Minimal k-mer length for assembly
#' @param maxKmer Maximal k-mer length for assembly. 0 means auto
#' @param steps Number of assembly iterations from minimal to maximal k-mer length
#' @param minCount Minimal count for k-mers retained
#' @param maxKmerCount Maximum k-mer count for fork tie-breaking
#' @param vectorPercent Percentage of reads containing 19-mer for adapter detection. 1.0 disables
#' @param insertSize Expected insert size for paired reads. 0 means auto
#' @param fraction Maximum noise to signal ratio acceptable for extension
#' @param maxSnpLen Maximal SNP length
#' @param minContig Minimal contig length reported in output
#' @param allowSnps Allow additional step for SNP discovery
#' @param forceSingleEnds Do not use paired-end information
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#'
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetMapCellSKESAintegrated <- function(
    bascetRoot,
    inputName="filtered",
    outputName="contigs",
    numThreads=NULL,
    numSkesaWorkers=NULL,
    numSkesaCores=NULL,
    numThreadsRead=NULL,
    totalMem=NULL,
    kmer=21,
    maxKmer=0,
    steps=11,
    minCount=1,
    maxKmerCount=10,
    vectorPercent=0.05,
    insertSize=0,
    fraction=0.01,
    maxSnpLen=150,
    minContig=50,
    allowSnps=FALSE,
    forceSingleEnds=FALSE,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  if(is.null(numThreads)) {
    numThreads <- as.integer(runner@ncpu)
  }

  if(is.null(numSkesaWorkers) && is.null(numSkesaCores)) {
    numSkesaWorkers <- 1
    numSkesaCores <- numThreads
  } else if(is.null(numSkesaWorkers)) {
    numSkesaWorkers <- max(1, floor(numThreads / numSkesaCores))
  } else if(is.null(numSkesaCores)) {
    numSkesaCores <- max(1, floor(numThreads / numSkesaWorkers))
  }

  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.valid.threadcount(numThreads))
  stopifnot(is.valid.threadcount(numSkesaWorkers))
  stopifnot(is.valid.threadcount(numSkesaCores))
  if(!is.null(numThreadsRead)) {
    stopifnot(is.valid.threadcount(numThreadsRead))
  }
  stopifnot(is.positive.integer(kmer))
  stopifnot(is.integer.like(maxKmer), maxKmer >= 0)
  stopifnot(is.positive.integer(steps))
  if(!is.null(minCount)) {
    stopifnot(is.positive.integer(minCount))
  }
  stopifnot(is.positive.integer(maxKmerCount))
  stopifnot(is.numeric(vectorPercent), vectorPercent >= 0, vectorPercent <= 1)
  stopifnot(is.integer.like(insertSize), insertSize >= 0)
  stopifnot(is.numeric(fraction), fraction > 0)
  stopifnot(is.positive.integer(maxSnpLen))
  stopifnot(is.positive.integer(minContig))
  stopifnot(is.logical(allowSnps))
  stopifnot(is.logical(forceSingleEnds))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))

  #Check memory sizes
  totalMem <- checkTotalMemArg(totalMem, runner, bascetInstance)

  inputFiles <- file.path(bascetRoot, detectShardsForFile(bascetRoot, inputName))
  num_shards <- length(inputFiles)

  if(num_shards==0){
    stop("No input files")
  }

  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "zip", num_shards)

  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    RunJob(
      runner = runner,
      jobname = "Zskesa",
      bascetInstance = bascetInstance,
      cmd = c(
        shellscriptMakeBashArray("files_in",inputFiles),
        shellscriptMakeBashArray("files_out",outputFiles),

        ### Abort early if needed
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),

        assembleBascetCommand(bascetInstance, c(
          "skesa",
          "-i=${files_in[$TASK_ID]}",
          "-o=${files_out[$TASK_ID]}",
          paste0("--skesa-workers=",numSkesaWorkers),
          paste0("--skesa-cores=",numSkesaCores),
          if(!is.null(numThreadsRead)) paste0("--num-threads-read=",numThreadsRead),
          if(!is.null(totalMem)) paste0("--memory=",format_size_bascet(totalMem)),
          paste0("--kmer=",kmer),
          paste0("--max-kmer=",maxKmer),
          paste0("--steps=",steps),
          if(!is.null(minCount)) paste0("--min-count=",minCount),
          paste0("--max-kmer-count=",maxKmerCount),
          paste0("--vector-percent=",vectorPercent),
          paste0("--insert-size=",insertSize),
          paste0("--fraction=",fraction),
          paste0("--max-snp-len=",maxSnpLen),
          paste0("--min-contig=",minContig),
          if(allowSnps) "--allow-snps",
          if(forceSingleEnds) "--force-single-ends",
          "--single-pass-counter"
        ))
      ),
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
}


################################################################################
################ QUAST #########################################################
################################################################################



###############################################
#' Run QUAST on reads of all cells.
#' This is a thin wrapper around BascetMapCell
#'
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param ... Additional arguments passed to \code{\link{BascetMapCell}}
#'
#' @return A job to be executed, or being executed, depending on runner settings
#' @seealso \code{\link{BascetMapCell}}
#' @export
BascetMapCellQUAST <- function(
    bascetRoot,
    inputName="filtered",
    outputName="quast",
    ...
){
  BascetMapCell(
    bascetRoot=bascetRoot,
    withfunction="_quast",
    inputName=inputName,
    outputName=outputName,
    ...
  )
}

###############################################
#' Callback function for aggregating QUAST data.
#' To be called from BascetAggregateMap
#' 
#' @param bascetFile An opened Bascet file
#' @param cellID Cell ID
#' @param bascetInstance A Bascet instance
#' 
#' @return Data for one cell
#' @export
aggr.quast <- function(
    bascetFile, 
    cellID, 
    bascetInstance
){
  
  #print(cellID)
  fcont <- BascetReadFile(
    bascetFile, 
    cellID, 
    "transposed_report.tsv", 
    as="text", 
    bascetInstance=bascetInstance
  )
  if(!is.null(fcont)){
    dat <- data.frame(
      row.names=stringr::str_split(fcont[1],"\t")[[1]],
      value=stringr::str_split(fcont[2],"\t")[[1]]
    )
    dat <- dat[-1,,drop=FALSE]
    
    rownames(dat) <- stringr::str_replace_all(rownames(dat), stringr::fixed("#"),"Number of")
    
    #Arrange in the right format
    dat <- t(dat)
    
    #TODO should set data types to double whenever possible
    
    dat    
  } else {
    data.frame()
  }
}



###############################################
#' Aggregate data from QUAST
#' This is a thin wrapper around BascetAggregateMap
#'
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param includeCells Character vector of cell names to include, or NULL for all cells
#' @param ... Additional arguments passed to \code{\link{BascetAggregateMap}}
#'
#' @return Aggregated data
#' @seealso \code{\link{BascetAggregateMap}}
#' @export
BascetAggregateQUAST <- function(
    bascetRoot,
    inputName="quast",
    includeCells=NULL,
    ...
){
  BascetAggregateMap(
    bascetRoot,
    inputName,
    aggr.quast,
    includeCells=includeCells,
    ...
  )
}






################################################################################
################ FASTQC ########################################################
################################################################################


###############################################
#' Run FASTQC on reads of all cells.
#'
#' This uses the bascet integrated fastqc command rather than the old mapcell
#' system.
#'
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param numThreads Total thread budget. Defaults to the runner CPU count
#' @param numThreadsRead Threads used by the TIRP reader. If NULL, use the CLI default
#' @param nogroup Do not group bases in the FastQC per-base modules
#' @param expgroup Use exponential base grouping in the FastQC per-base modules
#' @param kmerSize K-mer size for FastQC k-mer content
#' @param minLength Minimum sequence length to include
#' @param dupLength Length to truncate sequences for duplication detection
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#'
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetMapCellFASTQC <- function(
    bascetRoot,
    inputName="filtered",
    outputName="fastqc",
    numThreads=NULL,
    numThreadsRead=NULL,
    nogroup=FALSE,
    expgroup=FALSE,
    kmerSize=7,
    minLength=0,
    dupLength=50,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  if(is.null(numThreads)) {
    numThreads <- as.integer(runner@ncpu)
  }

  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.valid.threadcount(numThreads))
  if(!is.null(numThreadsRead)) {
    stopifnot(is.valid.threadcount(numThreadsRead))
    stopifnot(numThreads > numThreadsRead)
  }
  stopifnot(is.logical(nogroup))
  stopifnot(is.logical(expgroup))
  stopifnot(is.positive.integer(kmerSize))
  stopifnot(is.integer.like(minLength), minLength >= 0)
  stopifnot(is.positive.integer(dupLength))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))

  inputFiles <- file.path(bascetRoot, detectShardsForFile(bascetRoot, inputName))
  num_shards <- length(inputFiles)

  if(num_shards==0){
    stop("No input files")
  }

  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "zip", num_shards)

  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    RunJob(
      runner = runner,
      jobname = "Zfastqc",
      bascetInstance = bascetInstance,
      cmd = c(
        shellscriptMakeBashArray("files_in",inputFiles),
        shellscriptMakeBashArray("files_out",outputFiles),

        ### Abort early if needed
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),

        assembleBascetCommand(bascetInstance, c(
          "fastqc",
          "-i=${files_in[$TASK_ID]}",
          "-o=${files_out[$TASK_ID]}",
          paste0("--threads=",numThreads),
          if(!is.null(numThreadsRead)) paste0("--num-threads-read=",numThreadsRead),
          if(nogroup) "--nogroup",
          if(expgroup) "--expgroup",
          paste0("--kmer-size=",kmerSize),
          paste0("--min-length=",minLength),
          paste0("--dup-length=",dupLength)
        ))
      ),
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
}




###############################################
#' Aggregate data from FASTQC
#' This is a thin wrapper around BascetAggregateMap
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param includeCells Character vector of cell names to include, or NULL for all cells
#' @param ... Additional arguments passed to \code{\link{BascetAggregateMap}}
#'
#' @return Aggregated data
#' @seealso \code{\link{BascetAggregateMap}}
#' @export
BascetAggregateFASTQC <- function(
    bascetRoot,
    inputName="fastqc",
    includeCells=NULL,
    ...
){
  BascetAggregateMap(
    bascetRoot,
    inputName,
    aggr.fastqc,
    includeCells=includeCells,
    ...
  )
  #print(666)
  #CountDataFrameToSparseMatrix(m)
}



###############################################
#' Callback function for aggregating FASTQC data for each cell.
#' To be called from BascetAggregateMap
#' 
#' @param bascetFile An opened Bascet file
#' @param cellID Cell ID
#' @param bascetInstance A Bascet instance
#' 
#' @return Data for one cell
#' @export
aggr.fastqc <- function(
    bascetFile, 
    cellID, 
    bascetInstance
){
  
  #tmp <- readLines("/home/mahogny/github/zorn/test_aggr/abricate/stupid.tsv")
  #tmp <- readLines("/home/mahogny/github/zorn/test_aggr/abricate/salmonella_SRR33219394_ncbi.tsv")
  #lines <- readLines("/home/mahogny/github/zorn/test_aggr/fastqc/new/r1_fastqc_data.txt")
  #internal_parse_fastqc_data(lines)
  
  #Handle both R1 and R2
  lines_r1 <- BascetReadFile(bascetFile, cellID, "r1_fastqc_data.txt", as="text", bascetInstance=bascetInstance)
  lines_r2 <- BascetReadFile(bascetFile, cellID, "r2_fastqc_data.txt", as="text", bascetInstance=bascetInstance)

  
  #lines_r1 <- readLines("/home/mahogny/github/zorn/test_aggr/fastqc/new/r1_fastqc_data.txt")

  dat_r1 <- internal_parse_fastqc_data(lines_r1)
  dat_r2 <- internal_parse_fastqc_data(lines_r2)

  if(!is.null(dat_r1) | !is.null(dat_r2)){    
    #Return both R1 and R2
    list(
      r1=dat_r1, 
      r2=dat_r2
    )
  } else {
    print(paste("fastqc lengths", length(lines_r1), length(lines_r2)))
    NULL
  }
}


#Parse one FASTQC report file, divide into tables for each section
internal_parse_fastqc_data <- function(lines){
  
  if(length(lines)>0){
    #Separate state of section
    all_section_states <- list()
    
    #Put each separate section in a list
    list_sections <- list()
    lines <- lines[lines!=">>END_MODULE" & lines!=""]
    section_start <- which(stringr::str_detect(lines,stringr::fixed(">>")))
    section_end <- c(section_start[-1], length(lines))
    for(i in seq_along(section_start)) {
      subsection <- lines[section_start[i]:(section_end[i]-1)]

      if(length(subsection)>1){
        dat <- read.delim(text = subsection[-1])
      } else {
        dat <- data.frame()
      }
      
      #Divide section and state
      section_name_state <- stringr::str_sub(subsection[1],3)
      section_name_state <- stringr::str_split_fixed(section_name_state, "\t",2)
      
      all_section_states[[i]] <- data.frame(state=section_name_state[,1], val=section_name_state[,2])

      #Store this data frame in current section
      section_name <- section_name_state[1]
      list_sections[[section_name]] <- dat
    }
    
    all_section_states <- do.call(rbind,all_section_states)
    all_section_states <- data.frame(
      row.names=all_section_states$state,
      val=all_section_states$val
    )
    
    
    all_section_states <- t(all_section_states)
    #print(all_section_states)
    
    #Name section-state table and include it as well  
    #colnames(all_section_states) <- c("section","state")
    #all_section_states <- as.data.frame(all_section_states)
    list_sections[["section_state"]] <- all_section_states
    
    # list_sections <- internal_parse_fastqc_data(readLines("/home/mahogny/github/zorn/test_aggr/fastqc/new/r1_fastqc_data.txt"))
    
    #print(666)
    #print(all_section_states)
    
    list_sections    
  } else {
    NULL
  }
}




###############################################
#' Show the FASTQC HTML report for a cell, in the available web browser
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param cellID Name of the cell
#' @param readnum 1 or 2, for R1 and R2
#' @param useBrowser Use operating system browser to open file
#' @param verbose Show debug output
#' 
#' @return Nothing
#' @export
ShowFASTQCforCell <- function(
    bascetRoot,
    inputName,
    cellID, 
    readnum=c(1,2),
    useBrowser=FALSE,
    verbose=FALSE
){
  #check arguments
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.character(cellID))
  readnum <- match.arg(readnum)
  stopifnot(is.logical(useBrowser))
  stopifnot(is.logical(verbose))

  #Open the bascet file, get the HTML report
  if(verbose){
    print("Creating extract streamer session")
  }
  bascetFile <- OpenBascet(bascetRoot, inputName)
  if(verbose){
    print("Extract streamer session ok")
  }
  lines <- BascetReadFile(bascetFile, cellID, paste0("r",readnum,"_fastqc.html"), as="text", bascetInstance=bascetInstance)
  #lines <- readLines("/home/mahogny/github/zorn/test_aggr/fastqc/new/r1_fastqc.html")
  CloseBascet(bascetFile)
  
  #RStudio is picky. File to show really has to be created this way
  dir <- tempfile()
  dir.create(dir)
  htmlFile <- file.path(dir, "index.html")
  writeLines(lines, htmlFile)

  #Show the file
  viewer <- getOption("viewer")
  if (!is.null(viewer) & !useBrowser) {
    viewer(htmlFile)
  } else {
    utils::browseURL(htmlFile)
  }
  Invisible()
}


###############################################  
# Helper for FASTQC: get a data.frame across all cells for a section and read ############### what happen??
###############################################
#' Get a data frame for one type of FASTQ statistics across across all cells
#'
#' @param aggrFastqData Aggregated FASTQC data list from BascetAggregateFASTQC
#' @param section Name of the FASTQC section to extract (string)
#' @param readnum 1 or 2, for R1 and R2
#'
#' @return TODO
#' @export
GetFASTQCassembledDF <- function(
    aggrFastqData,
    section,
    readnum=c(1,2)
) {
  #check arguments
  readnum <- match.arg(readnum)

  #helper function
  internal_fastqc_getread_in_list <- function(lst,readnum){
    lapply(lst, function(s) s[[paste0("r",readnum)]])
  }

  #helper function
  internal_fastqc_add_cellid_to_list <- function(lst){
    lapply(names(lst), function(x) {
      temp <- as.data.frame(lst[[x]]) #to be on the safe side; fixed one bug
      temp$cellID <- x
      temp
    })
  }
  
  list_oneread <- internal_fastqc_getread_in_list(aggrFastqData,readnum) #Get data for given read
  #print(list_oneread)
  
  list_oneread_section <- lapply(list_oneread, function(s) s[[section]]) #Get sections for each cell
  
  df_section <- do.call(rbind, internal_fastqc_add_cellid_to_list(list_oneread_section))
  df_section
}


###############################################
#' From aggregated FASTQC data, plot adapter content
#' 
#' @param aggrFastqData Aggregated FASTQ data
#' @param readnum 1 or 2, for R1 or R2
#' 
#' @return A ggplot object
#' @export
PlotFASTQCadapterContent <- function(
    aggrFastqData,
    readnum=c(1,2)
) {
  #check arguments
  readnum <- match.arg(readnum)

  
  df_section <- GetFASTQCassembledDF(aggrFastqData,"Adapter Content",readnum)
  xlab_unit <- unique(df_section$X.Position)
  
  df_section$any_adapter <- rowSums(df_section[,!(colnames(df_section) %in% c("X.Position","cellID"))])
  df_section$X.Position <- factor(df_section$X.Position, levels=xlab_unit) 
  
  #Possible subsets; but would need one plot window for each
  #Illumina.Universal.Adapter Illumina.Small.RNA.3..Adapter Illumina.Small.RNA.5..Adapter Nextera.Transposase.Sequence SOLID.Small.RNA.Adapter cellID
  
  ggplot(df_section, aes(X.Position, any_adapter, group=cellID, color=cellID)) + 
    geom_line() + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylab("% Adapters")  
}






###############################################
#' From aggregated FASTQC data, get basic statistics for overlay on UMAP etc
#' 
#' @param aggrFastqData Aggregated FASTQ data
#' @param readnum 1 or 2, for R1 or R2
#' 
#' @return A data.frame
#' @export
GetFASTQCbasicStats <- function(
    aggrFastqData, 
    readnum=c(1,2)
) {
  #check arguments
  readnum <- match.arg(readnum)

  df <- GetFASTQCassembledDF(aggrFastqData,"Basic Statistics",readnum)
  df <- df[df$X.Measure %in% c("Sequence length","%GC","Sequences flagged as poor quality"),]
  df.mat <- as.data.frame(reshape2::acast(df, cellID ~ X.Measure, value.var = "Value"))
  
  seqlen <- stringr::str_split_fixed(df.mat$`Sequence length`,"-",2)
  
  data.frame(
    row.names = rownames(df.mat),
    gc=as.numeric(df.mat$`%GC`),
    num_seq_poor_quality=as.numeric(df.mat$`Sequences flagged as poor quality`),
    seqlen_from=as.numeric(seqlen[,1]),
    seqlen_to=as.numeric(seqlen[,2])
  )
}


###############################################
#' From aggregated FASTQC data, get overall pass-fail statistics for overlay on UMAP etc
#' 
#' @param aggrFastqData Aggregated FASTQC data list from BascetAggregateFASTQC
#' @param readnum 1 or 2, for R1 or R2
#'
#' @return A data.frame
#' @export
GetFASTQCpassfailStats <- function(
    aggrFastqData,
    readnum=c(1,2)
) {
  #check arguments
  readnum <- match.arg(readnum)
  
  df <- GetFASTQCassembledDF(aggrFastqData,"section_state",readnum)
  rownames(df) <- df$cellID
  df <- df[colnames(df)!="cellID",drop=FALSE]  
  #Could do both r1 and r2
  df  
}





if(FALSE){
  list_sections <- internal_parse_fastqc_data(readLines("/home/mahogny/github/zorn/test_aggr/fastqc/new/r1_fastqc_data.txt"))
  
  testpair <- list(
    r1=list_sections,
    r2=list_sections
  )
  
  testlist <- list(
    "A"=testpair,
    "B"=testpair
  )
  
  # [1] "Basic Statistics"
  # [1] "Per base sequence quality"
  # [1] "Per sequence quality scores"
  # [1] "Per base sequence content"
  # [1] "Per sequence GC content"
  # [1] "Per base N content"
  # [1] "Sequence Length Distribution"
  # [1] "Sequence Duplication Levels"
  # [1] "Overrepresented sequences"
  # [1] "Adapter Content"
  
  # PlotFASTQCadapterContent(testlist,1)
}


################################################################################
################ Abricate ######################################################
################################################################################



###############################################
#' Callback function for aggregating ABRicate data for each cell.
#' To be called from BascetAggregateMap
#' 
#' @param bascetFile An opened Bascet file
#' @param cellID Cell ID
#' @param bascetInstance A Bascet instance
#' 
#' @return Data for one cell
#' @export
aggr.abricate <- function(
    bascetFile, 
    cellID, 
    bascetInstance
){
  tmp <- BascetReadFile(bascetFile, cellID, "abricate.tsv", as="text", bascetInstance=bascetInstance)
  
  #tmp <- readLines("/home/mahogny/github/zorn/test_aggr/abricate/stupid.tsv")
  #tmp <- readLines("/home/mahogny/github/zorn/test_aggr/abricate/salmonella_SRR33219394_ncbi.tsv")
  dat <- read.delim(text = tmp)

  if(nrow(dat)>0){
    dat$cellID <- cellID  #needed? could make this general  ---- think no longer used!
  }
  
  #dat <- stringr::str_split_fixed(dat,"\t",15)
  #colnames(dat) <- c("FILE","SEQUENCE","START","END","STRAND","GENE","COVERAGE","COVERAGE_MAP","GAPS","PERC_COVERAGE","PERC_IDENTITY","DATABASE","ACCESSION","PRODUCT","RESISTANCE")
  #dat <- as.data.frame(dat)
  #dat$cellID <- c("A","A","A","B","B","B","B")
  dat
}


###############################################
#' Run Abricate on contigs of all cells.
#' This is a thin wrapper around BascetMapCell
#'
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param db Name of ABRicate database to use (string)
#' @param args Additional arguments passed to Bascet (list)
#' @param overwrite Whether to overwrite existing output (logical)
#' @param runner A runner object (e.g. LocalRunner, SlurmRunner)
#' @param bascetInstance A Bascet instance
#'
#' @return A job to be executed, or being executed, depending on runner settings
#' @seealso \code{\link{BascetMapCell}}
#' @export
BascetMapCellAbricate <- function(
    bascetRoot,
    inputName="contigs",
    outputName="abricate",
    db="ncbi",
    args=list(),
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  BascetMapCell(
    bascetRoot=bascetRoot,
    withfunction="_abricate",
    inputName=inputName,
    outputName=outputName,
    args = c(list(DATABASE_DIR=db), args),
    overwrite=overwrite,
    runner=runner,
    bascetInstance=bascetInstance
  )
}



###############################################
#' List installed databases available for Abricate
#'
#' @param dbdir Path to database directory (string)
#' @param bascetInstance A Bascet instance
#'
#' @return List of database names
#' @export
ListDatabaseAbricate <- function(
    dbdir,
    bascetInstance=GetDefaultBascetInstance()
) {
  #check arguments
  stopifnot(is.bascet.instance(bascetInstance))
  
  ret <- system(
    paste(
      bascetInstance@prependCmd,
      "abricate --list"
    ),
    intern = TRUE
  )
  ret
}




###############################################
#' Aggregate data from Abricate
#' This is a thin wrapper around BascetAggregateMap
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param includeCells Character vector of cell names to include, or NULL for all cells
#' @param ... Additional arguments passed to \code{\link{BascetAggregateMap}}
#'
#' @return Aggregated data
#' @seealso \code{\link{BascetAggregateMap}}
#' @export
BascetAggregateAbricate <- function(
    bascetRoot,
    inputName="abricate",
    includeCells=NULL,
    ...
){
  CountDataFrameToSparseMatrix(MapCellMultiListAsDataFrame(BascetAggregateMap(
    bascetRoot,
    inputName,
    aggr.abricate,
    includeCells=includeCells,
    ...
  )), "cellID","GENE")
}


################################################################################
################ Bakta #########################################################
################################################################################


###############################################
#' Run Bakta on contigs of all cells.
#' This is a thin wrapper around BascetMapCell
#'
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param db Path to database
#' @param args Additional arguments passed to Bascet (list)
#' @param overwrite Whether to overwrite existing output (logical)
#' @param runner A runner object (e.g. LocalRunner, SlurmRunner)
#' @param bascetInstance A Bascet instance
#'
#' @return A job to be executed, or being executed, depending on runner settings
#' @seealso \code{\link{BascetMapCell}}
#' @export
BascetMapCellBakta <- function(
    bascetRoot,
    inputName="contigs",
    outputName="bakta",
    db,
    args=list(),
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  BascetMapCell(
    bascetRoot=bascetRoot,
    withfunction="_bakta",
    inputName=inputName,
    outputName=outputName,
    args = c(list(DATABASE_DIR=db), args),
    overwrite=overwrite,
    runner=runner,
    bascetInstance=bascetInstance
  )
}




###############################################
#' Download a database for Bakta
#'
#' @param dbdir Directory path to download database to (string)
#' @param dbtype Database type: "light" or "full"
#' @param bascetInstance A Bascet instance
#'
#' @export
DownloadDatabaseBakta <- function(
    dbdir,
    dbtype=c("light","full"),  #todo look up how to handle documentation for this
    bascetInstance=GetDefaultBascetInstance()
) {
  if(file.exists(dbdir)) {
    print("Database already exists; skipping")
  } else {
    system(
      paste(
        bascetInstance@prependCmd,
        "bakta_db download --output",
        dbdir,
        "--type",dbtype
      )
    )
  }
}



################################################################################
################ Ariba #########################################################
################################################################################


###############################################
#' Run Ariba on reads of all cells.
#' This is a thin wrapper around BascetMapCell
#'
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param db Path to database
#' @param args Additional arguments passed to Bascet (list)
#' @param overwrite Whether to overwrite existing output (logical)
#' @param runner A runner object (e.g. LocalRunner, SlurmRunner)
#' @param bascetInstance A Bascet instance
#'
#' @return A job to be executed, or being executed, depending on runner settings
#' @seealso \code{\link{BascetMapCell}}
#' @export
BascetMapCellAriba <- function(
    bascetRoot,
    inputName="filtered",
    outputName="ariba",
    db,
    args=list(),
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  BascetMapCell(
    bascetRoot=bascetRoot,
    withfunction="_ariba",
    inputName=inputName,
    outputName=outputName,
    args = c(list(DATABASE_DIR=db), args),
    overwrite=overwrite,
    runner=runner,
    bascetInstance=bascetInstance
  )
}




###############################################
#' Download database for Ariba
#' 
#' @param dbdir Directory in which to download database
#' @param ref Which reference to download
#' @param bascetInstance A Bascet instance
#' 
#' @export
DownloadDatabaseAriba <- function(
    dbdir,
    ref=c("ncbi", "argannot", "card",  "megares", "plasmidfinder", "resfinder", "srst2_argannot", "vfdb_core", "vfdb_full", "virulencefinder"),
    bascetInstance=GetDefaultBascetInstance()
) {
  #check arguments
  ref <- match.arg(ref)
  stopifnot(is.bascet.instance(bascetInstance))

  if(file.exists(dbdir)) {
    print("Database already exists; skipping")
  } else {
    dir.create(dbdir, recursive=TRUE)

    #generate temp files
    tmp.fa <- tempfile(fileext = ".fa")
    tmp.tsv <- tempfile(fileext = ".tsv") 
    out.prepareref <- file.path(dbdir,"out.prepareref")  #out.ncbi.prepareref
    
    # ariba getref ncbi out.ncbi
    # ariba prepareref -f out.ncbi.fa -m out.ncbi.tsv out.ncbi.prepareref
    
    ### Get ref
    system(
      paste(
        bascetInstance@prependCmd,
        "ariba getref ",ref, tmp, 
        dbdir
      )
    )
    
    ### prepare the database
    system(
      paste(
        bascetInstance@prependCmd,
        "ariba prepareref -f ",tmp.fa,
        "-m", tmp.tsv, out.prepareref
        #        dbdir
      )
    )
    
    stop("todo implement")
  }
  
  #todo if above is heavy, make it a slurm job
}


###############################################
#' Callback function for aggregating ARIBA data for each cell.
#' To be called from BascetAggregateMap
#' 
#' @param bascetFile An opened Bascet file
#' @param cellID Cell ID
#' @param bascetInstance A Bascet instance
#' 
#' @return Data for one cell
#' @export
aggr.ariba <- function(bascetFile, cellID, bascetInstance){
  tmp <- BascetReadFile(bascetFile, cellID, "report.tsv", as="text", bascetInstance=bascetInstance)
  #tmp <- readLines("/home/mahogny/github/zorn/test_aggr/ariba/report.tsv")
  dat <- read.delim(text = tmp)
  dat
}



###############################################
#' Aggregate data from Ariba
#' This is a thin wrapper around BascetAggregateMap
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param includeCells Character vector of cell names to include, or NULL for all cells
#' @param ... Additional arguments passed to \code{\link{BascetAggregateMap}}
#'
#' @return Aggregated data
#' @seealso \code{\link{BascetAggregateMap}}
#' @export
BascetAggregateAriba <- function(
    bascetRoot,
    inputName="ariba",
    includeCells=NULL,
    ...
){
  CountDataFrameToSparseMatrix(MapCellMultiListAsDataFrame(BascetAggregateMap(
    bascetRoot=bascetRoot,
    inputName=inputName,
    aggr.ariba,
    includeCells=includeCells,
    ...
  )), "cellID","cluster")
}



################################################################################
################ AMRfinder #####################################################
################################################################################



###############################################
#' Run AMRfinder on contigs of all cells.
#' This is a thin wrapper around BascetMapCell
#'
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param db Path to database
#' @param args Additional arguments passed to Bascet (list)
#' @param overwrite Whether to overwrite existing output (logical)
#' @param runner A runner object (e.g. LocalRunner, SlurmRunner)
#' @param bascetInstance A Bascet instance
#'
#' @return A job to be executed, or being executed, depending on runner settings
#' @seealso \code{\link{BascetMapCell}}
#' @export
BascetMapCellAMRfinder <- function(
    bascetRoot,
    inputName="contigs",
    outputName="AMRfinder",
    db,
    args=list(),
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  BascetMapCell(
    bascetRoot=bascetRoot,
    withfunction="_amrfinder",
    inputName=inputName,
    outputName=outputName,
    args = c(list(DATABASE_DIR=db), args),
    overwrite=overwrite,
    runner=runner,
    bascetInstance=bascetInstance
  )
}



###############################################
#' Download a database for AMRfinder
#' 
#' @param dbdir Directory to download database to
#' @param bascetInstance A Bascet instance
#' 
#' @export
DownloadDatabaseAMRfinder <- function(
    dbdir,
    bascetInstance=GetDefaultBascetInstance()
) {
  #check arguments
  stopifnot(is.bascet.instance(bascetInstance))
  
  if(file.exists(dbdir)) {
    print("Database already exists; skipping")
  } else {
    system(
      paste(
        bascetInstance@prependCmd,
        "amrfinder_update -d",
        dbdir
      )
    )
  }
}





###############################################
#' Callback function for aggregating ABRicate data for each cell.
#' To be called from BascetAggregateMap
#' 
#' @param bascetFile An opened Bascet file
#' @param cellID Cell ID
#' @param bascetInstance A Bascet instance
#' 
#' @return Data for one cell
#' @export
aggr.amrfinder <- function(
    bascetFile, 
    cellID, 
    bascetInstance
){
  tmp <- BascetReadFile(bascetFile, cellID, "amrfinder.tsv", as="text", bascetInstance=bascetInstance)

  dat <- tryCatch(
   read.delim(text = tmp),
   error = function(e) {
     warning(paste0("aggr.amrfinder: failed to parse output for cell ", cellID, ": ", e$message))
     return(c())
   }
 )

  dat
}


###############################################
#' Aggregate data from AMRfinder
#' This is a thin wrapper around BascetAggregateMap
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param includeCells Character vector of cell names to include, or NULL for all cells
#' @param getColumn Column name from AMRfinder output to use as feature (string)
#' @param ... Additional arguments passed to \code{\link{BascetAggregateMap}}
#'
#' @return Aggregated data
#' @seealso \code{\link{BascetAggregateMap}}
#' @export
BascetAggregateAMRfinder <- function(
    bascetRoot,
    inputName="AMRfinder",
    includeCells=NULL,
    getColumn="Element.symbol",  #Element.name is an option, if one also want full name
    ...
){
  inp <- BascetAggregateMap(
    bascetRoot=bascetRoot,
    inputName=inputName,
    aggr.amrfinder,
    includeCells=includeCells,
    ...
  )
  MapCellMultiListAsDataFrame(inp)
  #CountDataFrameToSparseMatrix(MapCellMultiListAsDataFrame(inp), "cellID",getColumn)
}



################################################################################
################ GECCO #########################################################
################################################################################




###############################################
#' Run GECCO on contigs of all cells.
#'
#' This uses the bascet integrated gecco command rather than the old mapcell
#' system.
#'
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param numThreads Total thread budget. Defaults to the runner CPU count
#' @param dataDir GECCO data directory containing HMM, CRF model, and InterPro files
#' @param threshold Minimum probability for cluster membership
#' @param cds Minimum number of annotated CDS in a cluster
#' @param noMask Do not mask ambiguous nucleotides during gene prediction
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#'
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetMapCellGECCO <- function(
    bascetRoot,
    inputName="contigs",
    outputName="gecco",
    numThreads=NULL,
    dataDir=NULL,
    threshold=0.8,
    cds=3,
    noMask=FALSE,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  if(is.null(numThreads)) {
    numThreads <- as.integer(runner@ncpu)
  }

  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.valid.threadcount(numThreads))
  if(!is.null(dataDir)) {
    stopifnot(dir.exists(dataDir))
  }
  stopifnot(is.numeric(threshold), threshold >= 0, threshold <= 1)
  stopifnot(is.positive.integer(cds))
  stopifnot(is.logical(noMask))
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))

  inputFiles <- file.path(bascetRoot, detectShardsForFile(bascetRoot, inputName))
  num_shards <- length(inputFiles)

  if(num_shards==0){
    stop("No input files")
  }

  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "zip", num_shards)

  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    RunJob(
      runner = runner,
      jobname = "Zgecco",
      bascetInstance = bascetInstance,
      cmd = c(
        shellscriptMakeBashArray("files_in",inputFiles),
        shellscriptMakeBashArray("files_out",outputFiles),

        ### Abort early if needed
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),

        assembleBascetCommand(bascetInstance, c(
          "gecco",
          "-i=${files_in[$TASK_ID]}",
          "-o=${files_out[$TASK_ID]}",
          paste0("--threads=",numThreads),
          if(!is.null(dataDir)) paste0("--data-dir=",dataDir),
          paste0("--threshold=",threshold),
          paste0("--cds=",cds),
          if(noMask) "--no-mask"
        ))
      ),
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
}


###############################################
#' Callback function for aggregating GECCO data for each cell.
#' To be called from BascetAggregateMap
#' 
#' @param bascetFile An opened Bascet file
#' @param cellID Cell ID
#' @param bascetInstance A Bascet instance
#' 
#' @return Data for one cell
#' @export
aggr.gecco <- function(
    bascetFile, 
    cellID, 
    bascetInstance
){ 
  
  tmp <- BascetReadFile(bascetFile, cellID, "clusters.tsv", as="text", bascetInstance=bascetInstance)
  dat <- read.delim(text = tmp)
  dat
}



###############################################
#' Aggregate data from GECCO
#' This is a thin wrapper around BascetAggregateMap
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param includeCells Character vector of cell names to include, or NULL for all cells
#' @param ... Additional arguments passed to \code{\link{BascetAggregateMap}}
#'
#' @return Aggregated data
#' @seealso \code{\link{BascetAggregateMap}}
#' @export
BascetAggregateGECCO <- function(
    bascetRoot,
    inputName="gecco",
    includeCells=NULL,
    ...
){
  BascetAggregateMap(
    bascetRoot=bascetRoot,
    inputName=inputName,
    aggr.gecco,
    includeCells=includeCells,
    ...
  )
}
