
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
  mem="character",
  direct="logical",
  verbose="logical",
  deleteScript="logical",
  benchmark="logical"
)
) 



#' @export
setClass("SlurmJob", slots=list(
  pid="character", 
  cmd="character", 
  logfile="character",
  jobname="character",
  arraysize="numeric"
)
) 


###################################################### TODO this causes errors in devtools::document() ; why?
#Caused by error in `current@target`:
#  ! no applicable method for `@` applied to an object of class "environment"

#setMethod("show", "SlurmJob", function(object) {
  # cat(
  #   paste0(
  #     "PID:",object@pid, "  ",
  #     "array_size:",object@arraysize, "  ",
  #     "cmd:",object@cmd, "  ",
  #     "logfile:",object@logfile
  #   )
  # )
#})


###############################################
#' Create a runner that submits jobs to SLURM
#' 
#' @param settings Default settings to override; can be NULL
#' @param ncpu Number of cores requested (SLURM -c)
#' @param partition Which partition to run the job on (SLURM -p)
#' @param account Which account to run the job on (SLURM -A)
#' @param time The time the job is allowed to run, e.g. "0-72:00:00" (SLURM -t)
#' @param prepend Something to prepend to the command. TODO seems not used. present in instance instead!!
#' @param mem Amount of main memory to reserve (SLURM --mem)
#' @param deleteScript Delete job script after execution. Set to FALSE if you want to dissect it for debugging purposes
#' @param direct Run and get the result directly. FALSE implies asynchronous execution
#' @param benchmark Enable logging of final CPU and memory usage
#' @param verbose Enable additional debug output
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
    mem=NULL,
    direct=NULL,
    deleteScript=NULL,
    benchmark=NULL,
    verbose=NULL
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
      mem="",
      direct=TRUE,
      deleteScript=TRUE,
      benchmark=FALSE,
      verbose=FALSE
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

  
  if(!is.null(direct)){
    settings@direct <- direct
  }
  
  if(!is.null(deleteScript)){
    settings@deleteScript <- deleteScript
  }

  if(!is.null(benchmark)){
    settings@benchmark <- benchmark
  }
  
  if(!is.null(verbose)){
    settings@verbose <- verbose
  }
  
  settings
}




#' @export
setMethod(
  f = "RunJob",
  signature ="SlurmRunner",
  definition = function(runner, jobname, bascetInstance, cmd, arraysize) {
    
    print("Running job with slurm")
    
    ## Set up SBATCH settings
    scriptcontent <- c("#!/usr/bin/env bash")
    if(runner@account!=""){
      scriptcontent <- c(scriptcontent, paste("#SBATCH -A",runner@account))
    }
    if(runner@ncpu!=""){
      scriptcontent <- c(scriptcontent, paste("#SBATCH -c",runner@ncpu))  #number of cores https://docs.hpc2n.umu.se/documentation/batchsystem/batch_scripts/
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

    ## Decide on a tempdir location; different for each job
    scriptcontent <- c(
      scriptcontent,
      paste0("BASCET_TEMPDIR=",file.path(GetBascetTempDir(bascetInstance), "$SLURM_ARRAY_TASK_ID"))
    )
    
    ## Add the command
    this_cmd <- stringr::str_replace_all(cmd,stringr::fixed("$TASK_ID"),"$SLURM_ARRAY_TASK_ID")
    this_cmd <- stringr::str_replace_all(this_cmd,stringr::fixed("${TASK_ID}"),"${SLURM_ARRAY_TASK_ID}")
    scriptcontent <- c(scriptcontent, this_cmd)
    
    scriptcontent <- c(scriptcontent, "echo \"End of SLURM script\"")
    
    
    ## Write the script to a temporary file
    if(runner@deleteScript) {
      slurm_script <- tempfile(fileext = ".sh")
    } else {
      slurm_script <- "bascet_last_slurm_script.sh"
    }
    
    #print(slurm_script)    
    writeLines(con=slurm_script, scriptcontent)
    
    if(runner@verbose){
      #/var/spool/slurmd/job34516237/slurm_script: line 9: [: missing `]'   --- need to fix
      writeLines("=========== SLURM script start ===========")
      writeLines(scriptcontent)
      writeLines("=========== SLURM script end =============")
    }
        
    ## Run the script; catch the message from sbatch
    cmd <- paste0("sbatch --array=0-",arraysize-1,"  ",slurm_script)
    ret <- system(cmd, intern = TRUE)
    print(ret)

    if(length(ret)==0){
      print("Failed to start job")
      print(slurm_script)
      return(NULL)
    }
    
    ## Check if it worked, or if there was an error
    if(stringr::str_starts(ret, "Submitted batch job ")){
      
      ## Remove the temporary file. Worst case, done at the end of the R session, but better done earlier
      if(runner@deleteScript) {
        file.remove(slurm_script)
      }
      
      pid <- stringr::str_remove(ret, stringr::fixed("Submitted batch job "))
      
      #Create the job
      job <- new(
        "SlurmJob",
        pid=pid,
        cmd=cmd, 
        logfile="todo",
        jobname=jobname,
        arraysize=arraysize
      )
      
      if(runner@direct){
        #Wait for the job if direct mode
        WaitForJob(job)
        new_no_job()
      } else {
        #Return the job with PID set  
        job
      }
      
    } else {
      stop("Failed to start job (err2)")
      #print(slurm_script)
      #NULL
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
    
    cli::cli_progress_bar(
      total = job@arraysize,
      format = "{cli::pb_bar} {cli::pb_percent} | {cur_summary}"
    )
    
    while(TRUE) {
      info <- JobStatus(job)
      
      if(nrow(info)==0){
        #print("Waiting to start")  
        cur_summary <- "Waiting for SLURM job to start"
        cli::cli_progress_update(set = 0)
        Sys.sleep(5)
      } else {
        if(all(info$status=="COMPLETED")){
          break;
        } else {
          #print(info$status)
          
          #TODO: list only jobs starting with #_#" "  ; no . in it!   then no need for /2 below
          
          
          num_total <- job@arraysize
          num_running <- floor(sum(info$status=="RUNNING")/2)
          num_completed <- floor(sum(info$status=="COMPLETED")/2)  ### for some reason, these are reported twice.. ish
          num_failed <- floor(sum(info$status=="FAILED")/2)  ### for some reason, these are reported twice.. ish ?
          num_outofmem <- floor(sum(stringr::str_starts(info$status,"OUT_OF_ME"))/2)  ### for some reason, these are reported twice.. ish ?
          num_cancelled <- floor(sum(info$status=="CANCELLED")/2)  ### for some reason, these are reported twice.. ish ?
          

          cur_summary <- paste0(
            job@pid," ",
            job@jobname,"   ",
            "Total to run: ",num_total,"   ",
            "Completed: ",num_completed,"   ",
            "Running: ",num_running,"   ",
            "Failed: ", num_failed,"   ",
            "Cancelled: ", num_cancelled,"   ",
            "Out-of-mem: ", num_outofmem
          )
          cli::cli_progress_update(set = num_completed)

          Sys.sleep(5)
        }
      }
    }
      
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
    
#    ret <- system(paste("sacct -j",job@pid), intern = TRUE)
#    print(ret)
        
    
    ret <- system(paste("sacct -j",job@pid, "-o state"), intern = TRUE)
    ret <- ret[-c(1:2)]  #This is " state " and "-----"
    status <- stringr::str_remove_all(ret, " ")  
    data.frame(status=status) #, chard=1
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

