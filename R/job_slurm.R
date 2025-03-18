
################################################################################
################ Things for running SLURM ######################################
################################################################################


#' @export
setClass("SlurmRunner", slots=list(
  settings="character", 
  ncpu="character", 
  partition="character", 
  account="character", 
  time="character",
  prepend="character",
  mem="character"
)
) 



#' @export
setClass("SlurmJob", slots=list(
  pid="character", 
  cmd="character", 
  logfile="character",
  arraysize="numeric"
)
) 



###############################################
#' Create a runner that submits jobs to SLURM
#' 
#' @param settings Default settings to override; can be NULL
#' @param ncpu description
#' @param partition description
#' @param account description
#' @param time The time the job is allowed to run, e.g. "0-72:00:00"
#' @param prepend description
#' @param mem description
#' 
#' @return A SLURM runner
#' @export
SlurmRunner <- function(
    settings=NULL, 
    ncpu=NULL, 
    partition=NULL, 
    account=NULL, 
    time=NULL, 
    prepend=NULL, 
    mem=NULL
){
  
  ## Create a new default
  if(is.null(settings)){
    settings <- new(
      "SlurmRunner", 
      ncpu="1", 
      partition="",  #or NULL?? 
      account="",  #or NULL??
      time="0-72:00:00",
      prepend="",
      mem=""
    )
  }
  
  ##add new stuff
  if(!is.null(ncpu)){
    settings@ncpu <- ncpu
  }
  
  if(!is.null(partition)){
    settings@partition <- partition
  }
  
  if(!is.null(account)){
    settings@account <- account
  }
  
  if(!is.null(time)){
    settings@time <- time
  }
  
  if(!is.null(prepend)){
    settings@prepend <- prepend
  }
  
  if(!is.null(mem)){
    settings@mem <- mem
    #mem can have "g" at end. need to parse to integer  TODO
  }
  
  settings
}




#' @export
setMethod(
  f = "RunJob",
  signature ="SlurmRunner",
  definition = function(runner, jobname, cmd, arraysize) {
    
    print("Running job with slurm")
    
    ## Set up SBATCH settings
    scriptcontent <- c("#!/usr/bin/env bash")
    if(runner@account!=""){
      scriptcontent <- c(scriptcontent, paste("#SBATCH -A",runner@account))
    }
    if(runner@ncpu!=""){
      scriptcontent <- c(scriptcontent, paste("#SBATCH -n",runner@ncpu))
    }
    if(runner@time!=""){
      scriptcontent <- c(scriptcontent, paste("#SBATCH -t",runner@time))
    }
    if(runner@partition!=""){
      scriptcontent <- c(scriptcontent, paste("#SBATCH -p",runner@partition))
    }
    if(runner@mem!=""){
      scriptcontent <- c(scriptcontent, paste("#SBATCH --mem",runner@mem))
    }
    scriptcontent <- c(scriptcontent, paste("#SBATCH -J ",jobname))

    ## Add the command
    this_cmd <- stringr::str_replace_all(cmd,stringr::fixed("$TASK_ID"),"$SLURM_ARRAY_TASK_ID")
    scriptcontent <- c(scriptcontent, this_cmd)
    
    ## Write the script to a temporary file
    slurm_script <- tempfile(fileext = ".sh")
    print(slurm_script)    
    writeLines(con=slurm_script, c(
      "#!/bin/bash",
      scriptcontent)
    )
    print(scriptcontent)
    
    ## Run the script; catch the message from sbatch
    cmd <- paste0("sbatch --array=0-",arraysize-1,"  ",slurm_script)
    ret <- system(cmd, intern = TRUE)
    print(ret)
    
    ## Remove the temporary file. Worst case, done at the end of the R session, but better done earlier
    file.remove(slurm_script)
    
    
    ## Check if it worked, or if there was an error
    if(stringr::str_starts(ret, "Submitted batch job ")){
      pid <- stringr::str_remove(ret, stringr::fixed("Submitted batch job "))
      
      #Return the job with PID set  
      job <- new(
        "SlurmJob",
        pid=pid,
        cmd=cmd, 
        logfile="todo",
        arraysize=arraysize
      )
    } else {
      stop("Failed to start job")
      NULL
    }
  }
)





#Has possibility of ctrl+c; just keeps polling, possibly with a status indicator from log. or keep plotting log file
#' @export
setMethod(
  f = "WaitForJob",
  signature ="SlurmJob",
  definition = function(job) {
    # sacct -j JOBID -o jobid,submit,start,end,state
    while(TRUE) {
      info <- JobStatus(job)
      if(info$status=="RUNNING" || info$status=="PENDING"){
        print(info$status)
        Sys.sleep(5)
      } else {
        break;
      }
    }
    print(info$status)    
  }
)



#' @export
setMethod(
  f = "CancelJob",
  signature ="SlurmJob",
  definition = function(job) {
    #First check if it is running at all
    
    system(paste("scancel", job@pid))
  }
)




#' @export
setMethod(
  f = "JobStatus",
  signature ="SlurmJob",
  definition = function(job) {
    # sacct -j JOBID -o jobid,submit,start,end,state
    
    #ret <- system(paste("sacct -j",job@pid, "-o jobid,submit,start,end,state"), intern = TRUE)
    #print(ret)
    ret <- system(paste("sacct -j",job@pid, "-o state"), intern = TRUE)[3]
    status <- stringr::str_remove_all(ret, " ")  
    data.frame(status=status, chard=1)
    
    
    ## not enough! jobarray
  }
)


#' @export
setMethod(
  f = "JobLog",
  signature ="SlurmJob",
  definition = function(job) {
    #read this file
    
    job@logfile
  }
)



####### testing
if(FALSE){
  
  inst <- SlurmRunner(partition="shared", account = "naiss2024-22-1647", ncpu = "10")
  job <- RunJob(inst,"ls $JOBID",1)
  
  JobStatus(job)
  CancelJob(job)
  
  
  
  ret <- system(paste("sacct -j",job@pid, "-o jobid,submit,start,end,state"), intern = TRUE)
  ret <- system(paste("sacct -j",job@pid, "-o jobid,state"), intern = TRUE)[3]
  ret <- system(paste("sacct -j",job@pid, "-o state"), intern = TRUE)[3]
  status <- stringr::str_remove_all(ret, " ")  
  #first line after header
  ret[3]
}

