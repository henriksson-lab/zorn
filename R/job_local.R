

################################################################################
################ Things for running locally ####################################
################################################################################


###############################################
#' @export
setClass("LocalRunner", slots=list(
  ncpu="character",
  mem="character",
  direct="logical",
  showScript="logical"
)
) 


###############################################
#' @export
setClass("LocalJob", slots=list(
  cmd="character",
  proc="ANY", #processx:process
  logfile="character",
  arraysize="numeric"
)
) 

###############################################
#' Create new local runner instance
#' 
#' @param settings Default settings to override; can be NULL
#' @param ncpu Number of CPU cores to use. By default, will attempt to detect and use all of them
#' @param mem Total memory for each job. This setting will be used whenever possible. Default is to use up to all available memory
#' @param direct Run jobs synchronously
#' @param showScript Show the script to run, for debugging purposes
#' 
#' @return A runner instance
#' @export
LocalRunner <- function(
    settings=NULL,
    ncpu=NULL, 
    mem=NULL,
    direct=TRUE, 
    showScript=FALSE
){
  #Bring in defaults
  if(!is.null(settings) & is.null(cpu)) {
    ncpu <- settings@ncpu
  }

  if(!is.null(settings) & is.null(mem)) {
    mem <- settings@mem
  }

  if(!is.null(settings) & is.null(direct)) {
    direct <- settings@direct
  }
  
  if(!is.null(settings) & is.null(showScript)) {
    showScript <- settings@showScript
  }
  
    
  #Try to detect parameters    
  if(is.null(ncpu)) {
    ncpu <- as.character(parallel::detectCores())
    if(is.na(ncpu)) {
      warning("Cannot detect number of CPU cores so it was set to 1. Set it manually instead")
      ncpu <- "1"
    }
  }
  
  
  if(is.null(mem)) {
    bram <- benchmarkme::get_ram()
    if(is.null(bram)) {
      print("Warning: unable to detect amount of ram; it is better if you specify it. Setting a default of 64g")
      mem <- "64g"
    } else {
      bram_gb <- round(as.double(bram)/1000000000)
      mem <- paste0(bram_gb,"g")
    }
  } else {
    #Check that memory can be parsed and is some bare minimum
    stopifnot(parse_size_to_mb(mem)>1000)
  }
  
  new(
    "LocalRunner", 
    ncpu=ncpu, 
    mem=mem,
    direct=direct, 
    showScript=showScript
  )
}




###############################################
#' @export
setMethod(
  f = "RunJob",
  signature ="LocalRunner",
  definition = function(runner, jobname, bascetInstance, cmd, arraysize) { 
    print("Starting local job")
    
    ## Decide on a tempdir location; different for each job. Ensure to create it
    path_tempdir <- file.path(GetBascetTempDir(bascetInstance), "$TASK_ID")
    cmd <- c(
      paste0("BASCET_TEMPDIR=",path_tempdir),
      paste0("mkdir -p ",path_tempdir),
      cmd
    )

    ## Decide on a log location; different for each job. Create the directory
    cmd <- c(
      "mkdir -p logs",
      paste0("BASCET_LOGFILE=",paste0("logs/",jobname,".${TASK_ID}.log")),
      cmd
    )
    
    cmd <- stringr::str_flatten(cmd,"\n")
    
    print("----------- arg")
    if(runner@showScript){
      print("=============== final local script start =================")
      writeLines(cmd)
      print("=============== final local script end =================")
    }
    
    
    #### Concatenate all scripts over all TASK_ID
    all_cmd <- c()
    for(i in seq_len(arraysize)){
      this_cmd <- stringr::str_replace_all(cmd,stringr::fixed("$TASK_ID"),i-1)  #0... TASK_ID -1
      this_cmd <- stringr::str_replace_all(this_cmd,stringr::fixed("${TASK_ID}"),i-1)  #0... TASK_ID -1
      ##  onep <- processx::process$new(this_cmd)
      all_cmd <- c(all_cmd,this_cmd)
    }
    
    
    if(runner@showScript){
      print("=============== all script start =================")
      writeLines(all_cmd)
      print("=============== all script end =================")
    }
    
    #Figure out name of file to store in
    tfile <- tempfile(pattern = jobname, fileext=".sh") #putting jobname here helps it show up in "ps"; but may cause issues if bad jobname given
    
    #Delete the file upon exit
    all_cmd <- c(
      paste("trap \"rm -rf ",tfile,"\" EXIT"),
      all_cmd
    )
    
    
    #Write script file
    writeLines(con=tfile,c(
      "#!/bin/bash",
      stringr::str_flatten(all_cmd,"\n")
    ))
    
    if(runner@direct) {
      print("Running directly")  
      system(paste("bash", tfile))
      #print(tfile)
      #file.remove(tfile) #assumes process has started. can we do better? check https://www.linuxjournal.com/content/bash-trap-command 
      
      #Return no job, as it is done already
      new_no_job()
    } else {
      print("Using separate process")
      
      #write(tfile,paste("\n", "rm", tfile), append=TRUE)  ##this will make the script delete itself at the end ; or put command first? ---- can also 
      job <- new(
        "LocalJob",
        cmd=cmd,
        proc=processx::process$new("bash",tfile),
        logfile="foo",
        arraysize=arraysize
      )
      
      #Sys.sleep(1)
      #file.remove(tfile) #assumes process has started. can we do better? yes -- can put delete command into the tfile itself!! TODO
      
      #Return the job with PIDs set  
      job
    }
  }
)





###############################################
#' @export
setMethod(
  f = "CancelJob",
  signature ="LocalJob",
  definition = function(job) {
    if(job@proc$is_alive()){
      job@pid$kill()
    }
  }
)




###############################################
#' @export
setMethod(
  f = "JobStatus",
  signature ="LocalJob",
  definition = function(job) {
    if(job@proc$is_alive()){
      data.frame(status="RUNNING")
    } else {
      data.frame(status="DONE")  ### if killed, keep track and set cancelled
    }
  }
)


###############################################
#' @export
setMethod(
  f = "JobLog",
  signature ="LocalJob",
  definition = function(job) {
    #read this file
    #job@logfile
  }
)



################################################################################
########### Default runner #####################################################
################################################################################

#Should not set this! It will prevent GetDefaultBascetRunner() from
#searching outside the package
#bascetRunner.default <- LocalRunner(direct = TRUE, showScript=FALSE)



################################################################################
########### Testing ############################################################
################################################################################

if(FALSE){
  inst <- LocalRunner(direct=TRUE)
#  thejob <- RunJob(LocalRunner(),"ls",1)
  thejob <- RunJob(LocalRunner(),"ls",c(),"ls",1)
  
  CancelJob(thejob)
  JobStatus(thejob)
}