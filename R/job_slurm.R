
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
  benchmark="logical",
  logTime="logical"
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


######################################################
#this causes errors in devtools::document() ; why?
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
#' @param logTime Log execution time
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
    verbose=NULL,
    logTime=NULL
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
      verbose=FALSE,
      logTime=FALSE
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
    #Check that the argument is correct. will stop otherwise
    parse_size_to_bytes(mem)
    settings@mem <- mem
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
  
  if(!is.null(logTime)){
    settings@logTime <- logTime
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

    ## Get start time
    scriptcontent <- c(scriptcontent, "timeStart=$( date +\"%s\" )")


    ## Decide on a tempdir location; different for each job
    tmproot <- GetBascetTempDir(bascetInstance)
    scriptcontent <- c(
      scriptcontent,
      paste0("BASCET_TEMPDIR=",file.path(tmproot, ".", "$SLURM_ARRAY_TASK_ID"))
    )

    ## Also create the temp directory
    scriptcontent <- c(
      scriptcontent,
      paste("mkdir -p", file.path(tmproot, ".", "$SLURM_ARRAY_TASK_ID"))
    )
    
    ## Decide on a log location; different for each job
    scriptcontent <- c(
      scriptcontent,
      "mkdir -p logs",
      paste0("BASCET_LOGFILE=",paste0("logs/",jobname,".${SLURM_ARRAY_TASK_ID}.log"))
    )
    
    ## Add the provided command
    this_cmd <- stringr::str_replace_all(cmd,stringr::fixed("$TASK_ID"),"$SLURM_ARRAY_TASK_ID")
    this_cmd <- stringr::str_replace_all(this_cmd,stringr::fixed("${TASK_ID}"),"${SLURM_ARRAY_TASK_ID}")
    scriptcontent <- c(scriptcontent, this_cmd)

    ## Get end time
    scriptcontent <- c(scriptcontent, "timeEnd=$( date +\"%s\" )")
    
    ## Compute time spent running
    scriptcontent <- c(scriptcontent, "let timeDelta=$timeEnd-$timeStart")
    
    ## Log running time if requested
    if(runner@logTime) {
      scriptcontent <- c(scriptcontent, "echo \"$SLURM_JOB_NAME $SLURM_ARRAY_TASK_ID  $timeDelta\" >> bascet_timelog.txt")
    }
    
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




#' @export
setMethod(
  f = "JobStatus",
  signature ="SlurmJob",
  definition = function(job) {
    
    ret <- system(paste("sacct -j",job@pid, "-o jobid,state"), intern = TRUE)
    ret <- ret[-c(1:2)]  #This are the lines " jobid state " and "-----"
    
    #split columns
    info <- as.data.frame(stringr::str_split_fixed(stringr::str_trim(ret),"[ ]+",2))
    colnames(info) <- c("proc","status")
    
    #split jobid[..] and jobid_.. into two columns
    divproc <- stringr::str_split_fixed(info$proc,"[_\\.]",2)
    info$proc <- divproc[,1]
    info$num <- divproc[,2]
    
    #Correct cropped labels
    info$status[stringr::str_starts(info$status,"CANC")] <- "CANCELLED"
    info$status[stringr::str_starts(info$status,"RUN")] <- "RUNNING"
    info$status[stringr::str_starts(info$status,"FAIL")] <- "FAILED"
    info$status[stringr::str_starts(info$status,"OUT_OF_ME")] <- "OUT_OF_MEM"
    
    info
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
    
    #Because SLURM info is emptied overnight, we need to keep old info,
    #or it will get lost
    last_info <- NULL
    
    while(TRUE) {
      info <- JobStatus(job)[,c("status","num")]
      
      #Merge new info with last info. Only keep the latest entries
      info <- rbind(info, last_info)
      info <- info[!duplicated(info$num),] #remove older entries
      
      #For each expected job in the array, figure out it's status
      current_state <- rep("PENDING",job@arraysize) ## keep outside??
      
      if(nrow(info)==0) {
        #If we are extremely unlucky, this can suddenly happen overnight while waiting for jobs to be done

        #todo should likely keep status outside loop TODO
        
        #print("Waiting to start")  
        cur_summary <- "Waiting for SLURM job to start"
        cli::cli_progress_update(set = 0)
        Sys.sleep(5)
      } else {
        
        #Parse each line of info
        for(i in 1:nrow(info)) {
          this_num <- info$num[i]
#          print(info)
          if(this_num=="[0]") {
            #This is the entry for the batch array as a whole
            if(info$status=="CANCELLED") {
              current_state[1:job@arraysize] <- "CANCELLED"
            }
          } else {
            
            #check if string is numeric
            if(grepl("[0123456789]+$",this_num)) {
              
              #This is an individual index into the array
              slurm_index <- as.integer(this_num)
              current_state[slurm_index] <- info$status[i]
              
            } else {
              #this might be a weird entry such as 35556921_0.+  ; not sure what to make of these
            }
          }
        }        
        
        
        if(all(current_state=="COMPLETED")){
          break;
        } else {
          #print(info$status)
          
          #Count number of jobs of each type
          num_total <- nrow(info)
          num_running <- floor(sum(info$status=="RUNNING")/2)
          num_completed <- floor(sum(info$status=="COMPLETED")) 
          num_failed <- floor(sum(info$status=="FAILED"))
          num_outofmem <- floor(sum(info$status=="OUT_OF_MEM"))
          num_cancelled <- floor(sum(info$status=="CANCELLED"))
          
          #TODO: could also paste separate entries
          
          cur_summary <- paste0(
            job@pid," ",
            job@jobname,"   ",
            "To run: ",num_total,"   ",
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

# job <- createSlurmJobFromExisting(35556921,2)
# JobStatus(job)
# WaitForJob(job)




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
  f = "JobLog",
  signature ="SlurmJob",
  definition = function(job) {
    #read this file
    
    job@logfile
  }
)



###############################################
#' This creates a job object, linking to a running command.
#' Mainly used for development but can be used in case
#' a Zorn session died and you want to create a new
#' monitor
#' 
#' @param pid Process ID
#' @param arraysize Size of array (must be known)
#' @return A SLURM job object
#' 
#' @export
createSlurmJobFromExisting <- function(pid, arraysize) {
  new(
    "SlurmJob",
    pid=as.character(pid),
    cmd="prevcmd", 
    logfile="todo",
    jobname="prevjob",
    arraysize=arraysize
  )
}


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

