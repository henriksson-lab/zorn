
################################################################################
################ Things for running SLURM ######################################
################################################################################


setClass("SlurmInstance", slots=list(
  settings="character", 
  ncpu="character", 
  partition="character", 
  account="character", 
  time="character",
  prepend="character",
  mem="character"
)
) 



setClass("SlurmJob", slots=list(
  pid="character", 
  cmd="character", 
  logfile="character",
  arraysize="numeric"
)
) 



SlurmInstance <- function(settings=NULL, ncpu=NULL, partition=NULL, account=NULL, time=NULL, prepend=NULL, mem=NULL){
  
  ## Create a new default
  if(is.null(settings)){
    settings <- new(
      "SlurmInstance", 
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




setMethod(
  f = "RunJob",
  signature ="SlurmInstance",
  definition = function(runner, jobname, withdata, cmd, arraysize) {
    
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
    
    ## Add the data needed by the command
    scriptcontent <- c(scriptcontent, withdata)
    
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



setMethod(
  f = "CancelJob",
  signature ="SlurmJob",
  definition = function(job) {
    #First check if it is running at all
    
    system(paste("scancel", job@pid))
  }
)




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
  
  inst <- SlurmInstance(partition="shared", account = "naiss2024-22-1647", ncpu = "10")
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

