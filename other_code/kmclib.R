 


#' @export
setClass("Kmc3DB", slots=list(
  location="character"
)
) 



#' @export
KmcFromFile <- function(location){
  new(
    "Kmc3DB",
    location=location
  )
}




do.quote <- function(s) paste0("\"", s, "\"")




#' K-Mer Counter (KMC) ver. 3.2.1 (2022-01-04)
#' Usage:
#'   kmc [options] <input_file_name> <output_file_name> <working_directory>
#'   kmc [options] <@input_file_names> <output_file_name> <working_directory>
#'   Parameters:
#'   input_file_name - single file in specified (-f switch) format (gziped or not)
#' @input_file_names - file name with list of input files in specified (-f switch) format (gziped or not)
#' Options:
#'   -v - verbose mode (shows all parameter settings); default: false
#' -k<len> - k-mer length (k from 1 to 256; default: 25)
#' -m<size> - max amount of RAM in GB (from 1 to 1024); default: 12
#' -sm - use strict memory mode (memory limit from -m<n> switch will not be exceeded)
#' -hc - count homopolymer compressed k-mers (approximate and experimental)
#' -p<par> - signature length (5, 6, 7, 8, 9, 10, 11); default: 9
#' -f<a/q/m/bam/kmc> - input in FASTA format (-fa), FASTQ format (-fq), multi FASTA (-fm) or BAM (-fbam) or KMC(-fkmc); default: FASTQ
#' -ci<value> - exclude k-mers occurring less than <value> times (default: 2)
#' -cs<value> - maximal value of a counter (default: 255)
#' -cx<value> - exclude k-mers occurring more of than <value> times (default: 1e9)
#' -b - turn off transformation of k-mers into canonical form
#' -r - turn on RAM-only mode 
#' -n<value> - number of bins 
#' -t<value> - total number of threads (default: no. of CPU cores)
#' -sf<value> - number of FASTQ reading threads
#' -sp<value> - number of splitting threads
#' -sr<value> - number of threads for 2nd stage
#' -j<file_name> - file name with execution summary in JSON format
#' -w - without output
#' -o<kmc/kff> - output in KMC of KFF format; default: KMC
#' -hp - hide percentage progress (default: false)
#' -e - only estimate histogram of k-mers occurrences instead of exact k-mer counting
#' --opt-out-size - optimize output database size (may increase running time)
#' Example:
#'   kmc -k27 -m24 NA19238.fastq NA.res /data/kmc_tmp_dir/
#'   kmc -k27 -m24 @files.lst NA.res /data/kmc_tmp_dir/
#'   
#'   
#'   
#'   
#'   



### todo add the rest
#' @export
KmcBuild <- function(
    inputFile, 
    outputFile=tempfile(),   #but note... it makes new files which are not set as tempfiles! need to handle later
    workingDir=tempdir(), 
    kmerSize=25, 
    excludeBelow=2, 
    maxCounter=255, 
    excludeAbove=1000000000,
    outputFileFormat="kmc"    # -o<kmc/kff>
){
  
    #note: it can also take multiple input files! worth supporting later
  
    inputFile <- file.path(inputFile)
    bn <- basename(inputFile)
    
    #-f<a/q/m/bam/kmc> - input in FASTA format (-fa), FASTQ format (-fq), multi FASTA (-fm) or BAM (-fbam) or KMC(-fkmc); default: FASTQ
    file_format_flag <- ""
    if(stringr::str_detect(bn, stringr::fixed(".fa")) || stringr::str_detect(bn, stringr::fixed(".fasta"))) {
      file_format_flag <- "-fm" #this is multiple fasta. -fa is a single fasta. why two??
    } else if(stringr::str_detect(bn, stringr::fixed(".fq")) || stringr::str_detect(bn, stringr::fixed(".fastq"))) {
      file_format_flag <- "-fq"
    } else if(stringr::str_detect(bn, stringr::fixed(".bam"))) {
      file_format_flag <- "-fbam"
    } else {
      stop("Cannot tell file format")
    }
  
  
    delete_tempdir <- FALSE
    if(!file.exists(workingDir)){
      delete_tempdir <- TRUE
      dir.create(workingDir)
    }
  
 # cat(outputFile)
  
    cmd <- paste(
      "/home/mahogny/github/kmc/bin/kmc",
#      "kmc",
      file_format_flag,
      paste0("-k", kmerSize),
      paste0("-ci", excludeBelow),
      paste0("-cs", maxCounter),
      paste0("-o", outputFileFormat),  #-o<kmc/kff>
      do.quote(inputFile),
      do.quote(outputFile),
      do.quote(workingDir)
      ## above , format number properly -cx
      
      )
    #TODO check if it worked?

    cat(cmd)
    system(cmd)
    
    if(delete_tempdir){
      unlink(workingDir)
    }
    
    KmcFromFile(outputFile)
}




#' @export
KmcDump <- function(db, excludeBelow=NULL, excludeAbove=NULL){
  
  outp <- tempfile()
  cmd <- paste(
    "kmc_dump",
    db@location,
    outp,
    if(!is.null(excludeBelow)) paste0("-ci", excludeBelow),
    if(!is.null(excludeAbove)) paste0("-cx", excludeAbove)
  )
  #TODO check if it worked?
  
  #cat(cmd)
  system(cmd)
  
  #Fast parsing of file
  dat <- readLines(outp)
  dat <- stringr::str_split_fixed(dat, stringr::fixed("\t"),2)
  colnames(dat) <- c("kmer","count")
  dat <- as.data.frame(dat)
  dat$count <- as.integer(dat$count) 
  
  unlink(outp)
  
  dat
}


# 
# kmc_tools ver. 3.2.1 (2022-01-04)
# Usage:
#   kmc_tools [global parameters] <operation> [operation parameters]
# Available operations:
#   transform            - transforms single KMC's database
#   simple               - performs set operation on two KMC's databases
# complex              - performs set operation on multiple KMC's databases
#   filter               - filter out reads with too small number of k-mers
#  global parameters:
#   -t<value>            - total number of threads (default: no. of CPU cores)
#   -v                   - enable verbose mode (shows some information) (default: false)
#   -hp                  - hide percentage progress (default: false)
# Example:
# kmc_tools simple db1 -ci3 db2 -ci5 -cx300 union db1_union_db2 -ci10
# For detailed help of concrete operation type operation name without parameters:
# kmc_tools simple



# 
# Complex operation allows one to define operations for more than 2 input k-mers sets. Command-line syntax:
# kmc_tools complex <operations_definition_file>
#  operations_definition_file - path to file which define input sets and operations. It is text file with following syntax:
#  ______________________________________________________________________________ 
# |INPUT:                                                                        |
# |<input1>=<input1_db_path> [params]                                            |
# |<input2>=<input2_db_path> [params]                                            |
# |...                                                                           |
# |<inputN>=<inputN_db_path> [params]                                            |
# |OUTPUT:                                                                       |
# |<out_db>=<ref_input><oper[c_mode]><ref_input>[<oper[c_mode]><ref_input>[...]  |
# |[OUTPUT_PARAMS:                                                             __|
# |<output_params>]                                                           |  /
# |                                                                           | / 
# |___________________________________________________________________________|/  
# input1, input2, ..., inputN - names of inputs used to define equation
# input1_db, input2_db_path, ..., inputN_db_path - paths to k-mers sets
# For each input there are additional parameters which can be set:
#   -ci<value> - exclude k-mers occurring less than <value> times 
#   -cx<value> - exclude k-mers occurring more of than <value> times
# out_db_path - path to output database
# ref_input - one of input1, input2, ..., inputN
# oper - one of {*,-,~,+}, which refers to {intersect, kmers_subtract, counters_subtract, union}
# operator * has the highest priority. Other operators has equal priorities. Order of operations can be changed with parentheses
# for {*,~,+} it is possible to redefine counter calculation mode ([c_mode]). Available values: min, max, diff, sum, left, right (detailet description available in simple help message)
# output_params are:
#   -ci<value> - exclude k-mers occurring less than <value> times 
#   -cx<value> - exclude k-mers occurring more of than <value> times
#   -cs<value> - maximal value of a counter
# Example:
#  __________________________________________________________________ 
# |INPUT:                                                            |
# |set1 = kmc_o1 -ci5                                                |
# |set2 = kmc_o2                                                     |
# |set3 = kmc_o3 -ci10 -cx100                                      __|
# |OUTPUT:                                                        |  /
# |result = (set3 + min set1) * right set2                        | / 
# |_______________________________________________________________|/  
# 
# 


#### Could support formulas!

KmcOp <- function(
    list_db, 
    op, 
    outputFile=tempfile()){

  if(is.null(list_db)){
    stop("Input should be a named list")
  }
  
  ## Generate INPUT secrtion, containing list of input files
  script_content <- c("INPUT:")
  for(n in names(list_db)){
    script_content <- c(
      script_content, 
      paste(n, "=", list_db[[n]]@location)
      )
  }
  
  
  script_content <- c(script_content, "OUTPUT:")
  script_content <- c(script_content, paste0(outputFile,"=", op))
  
  ## Store script
  script_file <- tempfile()
  writeLines(script_content, con=script_file)  

  ## Run script
  cmd <- paste(
    "/home/mahogny/github/kmc/bin/kmc_tools",
#    "kmc_tools",
    "complex",
    script_file
  )
  #TODO check if it worked?
  print(cmd)
  system(cmd)

  ## Delete script
  #unlink(script_file)
  
  KmcFromFile(outputFile)
  
}









if(FALSE){

#  KmcBuild
    

  dbx <- KmcBuild("/husky/fromsequencer/240701_wgs_atcc1/trimmed/ref10/separate/Bacillus pacificus (ATCC 10987) CP086328.1.fasta", excludeBelow = 0, maxCounter = 3, outputFileFormat = "kff")  #this is allowed in kmc 3.2.4! but also .1??
  
  db1 <- KmcBuild("/husky/fromsequencer/240701_wgs_atcc1/trimmed/ref10/separate/Bacillus pacificus (ATCC 10987) CP086328.1.fasta", excludeBelow = 0, maxCounter = 1, outputFileFormat = "kff")  #this is allowed in kmc 3.2.4! but also .1??
  db2 <- KmcBuild("/husky/fromsequencer/240701_wgs_atcc1/trimmed/ref10/separate/Cereibacter sphaeroides (ATCC 17029) chr1.fasta", excludeBelow = 0, maxCounter = 1, outputFileFormat = "kff")
  db3 <- KmcBuild("/husky/fromsequencer/240701_wgs_atcc1/trimmed/ref10/separate/Cereibacter sphaeroides (ATCC 17029) chr1.fasta", excludeBelow = 0, maxCounter = 1, outputFileFormat = "kff")
  
  ## KFF and KMC take the same amount of space. KFF is a single file instead of two
  db2
  db3
  
  table(KmcDump(db1, excludeBelow = 0)$count)  #none with count <2

  ############ IMPORTANT!!! FLAGS MUST BE BEFORE FILE LIST!!!  
  
  list_db <- list()
  list_db[["bf"]] <- db1
  list_db[["cs"]] <- db2

  KmcOp(list_db, "bf+cs")
  merged_db <- KmcOp(list_db, "bf+cs")
  
  histo <- KmcDump(merged_db)
  head(histo)
  
  
  # kmc_tools transform /tmp/Rtmpmk0QGN/file3c088d20373178 set_counts 2  /tmp/out
  # ==> Error: kmc_tools currently does not support k-mer sets without counters. It will be implemented soon. If needed faster please contact authors or post an issue at https://github.com/refresh-bio/KMC/issues.

  #(base) mahogny@beagle:~$ kmc_tools transform /tmp/Rtmpmk0QGN/file3c088d20373178 compact  /tmp/out
  #Error: kmc_tools currently does not support k-mer sets without counters. It will be implemented soon. If needed faster please contact authors or post an issue at https://github.com/refresh-bio/KMC/issues.
  
  
  # API here; apt-get gives KMC 3.2.1
  # https://github.com/refresh-bio/KMC/wiki/Use-the-KMC-directly-from-code-through-the-API
  
  system(paste("kmc_tools transform", db1@location, "set_counts 2  /tmp/out"))  
  
  
  ################
  
  #{*,-,~,+}, which refers to {intersect, kmers_subtract, counters_subtract, union}
  
  path_all_genome <- "/husky/fromsequencer/240701_wgs_atcc1/trimmed/ref10/separate"
  list_db <- list()
  for(f in list.files(path_all_genome)){
    list_db[[f]] <- KmcBuild(file.path(path_all_genome, f), maxCounter = 1, excludeBelow = 0)
    
    #3.2.4 only:
    #kmc_core/kmc.h  ... update 11 months ago
    #Params.counter_max == 1
    #Warning: using counter_max == 1 will cause not storying counters in KMC output file, all counters will be assumed to be 1. This is experimental and is not currently supported in kmc_tools. Will be implemented soon.
    
    ## kmc_dump can read it; so merge-add should be simple?
    ## kmc_dump/kmc_dump.cpp

    
    ## set counters to 1??
    
    ## inline uint32 calc_counter_size(int64 cutoff_max, int64 counter_max)   here, returns 0
  }
  names(list_db) <- paste0("gen", 1:length(list_db))  
  sum_db <- KmcOp(list_db, stringr::str_flatten(names(list_db),"+"))
  sum_db <- KmcOp(list_db, stringr::str_flatten(names(list_db),"~"))
  
  
  #add up pairwise??
  
  #need way to binarize kmer count when checking for number; set presence count
  
  
  histo <- KmcDump(sum_db)
  table(histo$count)  # How unique are kmers among species?
  
  
  #TODO: meah. set max counter instead, then +
  
  #TODO: why is 2 so common? not 1?
    
  print(f)

    
  
  ##TODO: think we currently exclude below in bascet. this is a terrible idea!!!
  
}

