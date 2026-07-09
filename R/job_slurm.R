
################################################################################
################ Things for running SLURM ######################################
################################################################################


#' SLURM runner
#'
#' Runner backend that submits Bascet jobs to a SLURM scheduler.
#'
#' @slot settings Reserved settings string.
#' @slot ncpu Number of CPU cores requested.
#' @slot partition SLURM partition.
#' @slot account SLURM account.
#' @slot time SLURM time limit.
#' @slot prepend Command prefix.
#' @slot mem Memory limit string.
#' @slot direct Whether to wait for completion before returning.
#' @slot verbose Whether to print additional debug output.
#' @slot deleteScript Whether generated scripts are deleted after submission.
#' @slot benchmark Whether benchmark logging is enabled.
#' @slot logTime Whether task runtime logging is enabled.
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



#' SLURM job
#'
#' Job object for a submitted SLURM array job.
#'
#' @slot pid SLURM job identifier.
#' @slot cmd Submitted command.
#' @slot jobname Job name.
#' @slot arraysize Number of array tasks.
#' @export
setClass("SlurmJob", slots=list(
  pid="character", 
  cmd="character", 
  jobname="character",
  arraysize="numeric"
)
) 


slurm_normalize_state <- function(state) {
  state <- stringr::str_trim(state)
  state <- stringr::str_split_fixed(state, "[[:space:]]+", 2)[,1]

  state[stringr::str_starts(state, "CANC")] <- "CANCELLED"
  state[stringr::str_starts(state, "COMPLET")] <- "COMPLETED"
  state[stringr::str_starts(state, "CONFIGUR")] <- "CONFIGURING"
  state[stringr::str_starts(state, "FAIL")] <- "FAILED"
  state[stringr::str_starts(state, "OUT_OF_ME")] <- "OUT_OF_MEM"
  state[stringr::str_starts(state, "PEND")] <- "PENDING"
  state[stringr::str_starts(state, "RUN")] <- "RUNNING"
  state[stringr::str_starts(state, "TIMEOUT")] <- "TIMEOUT"

  state
}


slurm_terminal_states <- function() {
  c("COMPLETED", "FAILED", "OUT_OF_MEM", "CANCELLED", "TIMEOUT")
}


slurm_failure_states <- function() {
  setdiff(slurm_terminal_states(), "COMPLETED")
}


slurm_has_terminal_array_status <- function(info, arraysize) {
  if(nrow(info)==0) {
    return(FALSE)
  }
  if(arraysize <= 0) {
    return(FALSE)
  }

  whole_array_failed <- any(
    !is.na(info$num) &
      info$num=="[0]" &
      info$status %in% slurm_failure_states()
  )
  if(whole_array_failed) {
    return(TRUE)
  }

  info <- info[!is.na(info$num),,drop=FALSE]
  info <- info[grepl("^[0-9]+$", info$num),,drop=FALSE]
  if(nrow(info)==0) {
    return(FALSE)
  }

  slurm_index <- as.integer(info$num)
  info <- info[slurm_index >= 0 & slurm_index < arraysize,,drop=FALSE]
  if(nrow(info)==0) {
    return(FALSE)
  }

  expected <- as.character(seq.int(0, arraysize - 1))
  terminal <- unique(info$num[info$status %in% slurm_terminal_states()])
  all(expected %in% terminal)
}


slurm_expand_array_indices <- function(num) {
  num <- stringr::str_remove_all(num, "\\[|\\]")
  num <- stringr::str_remove(num, "%.*$")
  parts <- stringr::str_split(num, ",", simplify = FALSE)[[1]]
  indices <- character()

  for(part in parts) {
    part <- stringr::str_trim(part)
    if(grepl("^[0-9]+$", part)) {
      indices <- c(indices, part)
    } else if(grepl("^[0-9]+-[0-9]+$", part)) {
      bounds <- as.integer(stringr::str_split_fixed(part, "-", 2))
      indices <- c(indices, as.character(seq.int(bounds[1], bounds[2])))
    }
  }

  indices
}


slurm_parse_status_lines <- function(lines, separator = "|") {
  lines <- lines[nzchar(stringr::str_trim(lines))]
  if(length(lines)==0) {
    return(data.frame(proc=character(), status=character(), num=character()))
  }

  records <- list()
  for(line in lines) {
    fields <- stringr::str_split_fixed(stringr::str_trim(line), stringr::fixed(separator), 2)
    if(ncol(fields)<2 || fields[1,2]=="") {
      next
    }

    jobid <- fields[1,1]
    status <- slurm_normalize_state(fields[1,2])

    # Drop SLURM job steps. The array task state is represented by the
    # corresponding parent entry, while .batch/.extern rows can duplicate it.
    jobid <- stringr::str_remove(jobid, "\\..*$")
    jobid <- stringr::str_remove(jobid, "\\+$")

    parts <- stringr::str_split_fixed(jobid, "_", 2)
    proc <- parts[1,1]
    num <- parts[1,2]

    if(num=="") {
      records[[length(records)+1]] <- data.frame(
        proc=proc,
        status=status,
        num=NA_character_
      )
    } else if(grepl("^\\[", num)) {
      for(index in slurm_expand_array_indices(num)) {
        records[[length(records)+1]] <- data.frame(
          proc=proc,
          status=status,
          num=index
        )
      }
    } else if(grepl("^[0-9]+$", num)) {
      records[[length(records)+1]] <- data.frame(
        proc=proc,
        status=status,
        num=num
      )
    }
  }

  if(length(records)==0) {
    data.frame(proc=character(), status=character(), num=character())
  } else {
    do.call(rbind, records)
  }
}


slurm_parse_sacct_table <- function(lines) {
  lines <- lines[nzchar(stringr::str_trim(lines))]
  if(length(lines)<=2) {
    return(data.frame(proc=character(), status=character(), num=character()))
  }

  lines <- lines[-c(1:2)]
  records <- list()
  for(line in lines) {
    fields <- stringr::str_split_fixed(stringr::str_trim(line), "[[:space:]]+", 2)
    if(fields[1,1]=="" || fields[1,2]=="") {
      next
    }

    parsed <- slurm_parse_status_lines(
      paste0(fields[1,1], "|", fields[1,2]),
      separator = "|"
    )
    if(nrow(parsed)>0) {
      records[[length(records)+1]] <- parsed
    }
  }

  if(length(records)==0) {
    data.frame(proc=character(), status=character(), num=character())
  } else {
    do.call(rbind, records)
  }
}


slurm_system2 <- function(command, args) {
  tryCatch(
    system2(command, args = args, stdout = TRUE, stderr = TRUE),
    warning = function(w) character(),
    error = function(e) character()
  )
}


######################################################
#this causes errors in devtools::document() ; why?
#Caused by error in `current@target`:
#  ! no applicable method for `@` applied to an object of class "environment"

#setMethod("show", "SlurmJob", function(object) {
  # cat(
  #   paste0(
  #     "PID:",object@pid, "  ",
  #     "array_size:",object@arraysize, "  ",
  #     "cmd:",object@cmd
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
    #Check that the argument is correct
    stopifnot(is.valid.memsize(mem))
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




#' @describeIn RunJob SLURM runner method.
#' @export
setMethod(
  f = "RunJob",
  signature ="SlurmRunner",
  definition = function(runner, jobname, bascetInstance, cmd, arraysize) {
    
    cat("Running job with slurm\n")
    
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

    if(!is.jobscript(cmd)) {
      stop("SlurmRunner requires cmd to be a JobScript")
    }
    cmd <- renderJobScriptBash(cmd)
    
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
      #scriptcontent <- c(scriptcontent, paste0("let cpuName=\"\\\"", benchmarkme::get_cpu()$model_name, "\\\"\""))  ## check!!
      scriptcontent <- c(scriptcontent, "echo -e \"$SLURM_JOB_NAME\t$SLURM_ARRAY_TASK_ID\t$timeDelta\" >> bascet_timelog.txt")
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
      cat("Failed to start job!\n")
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




#' @describeIn JobStatus SLURM job method.
#' @export
setMethod(
  f = "JobStatus",
  signature ="SlurmJob",
  definition = function(job) {
    sacct_ret <- slurm_system2(
      "sacct",
      c("-j", job@pid, "--parsable2", "--noheader", "-o", "JobID,State")
    )
    sacct_info <- slurm_parse_status_lines(sacct_ret, separator = "|")

    if(nrow(sacct_info)==0) {
      sacct_ret <- slurm_system2(
        "sacct",
        c("-j", job@pid, "-o", "jobid,state")
      )
      sacct_info <- slurm_parse_sacct_table(sacct_ret)
    }

    if(slurm_has_terminal_array_status(sacct_info, job@arraysize)) {
      return(sacct_info)
    }

    # sacct can lag for just-submitted jobs or be unavailable on some clusters.
    # squeue sees pending/running tasks, including compact array ranges.
    squeue_ret <- slurm_system2(
      "squeue",
      c("-j", job@pid, "-h", "-o", "%i|%T")
    )
    squeue_info <- slurm_parse_status_lines(squeue_ret, separator = "|")

    info <- rbind(squeue_info, sacct_info)
    if(nrow(info)==0) {
      return(info)
    }

    info <- info[!is.na(info$num),,drop=FALSE]
    info$status <- slurm_normalize_state(info$status)
    info
  }
)





#Has possibility of ctrl+c; just keeps polling, possibly with a status indicator from log. or keep plotting log file
#' @describeIn WaitForJob SLURM job method.
#' @export
setMethod(
  f = "WaitForJob",
  signature ="SlurmJob",
  definition = function(job) {
    # sacct -j JOBID -o jobid,submit,start,end,state
    
    cur_summary <- "Waiting for SLURM job to start"
    cli::cli_progress_bar(
      total = job@arraysize,
      .auto_close = FALSE,
      format = "{cli::pb_bar} {cli::pb_percent} | {cur_summary}"
    )
    cli::cli_progress_update(set = 0)

    #Because SLURM info is emptied overnight, we need to keep old info, or it will get lost
    last_info <- NULL
    num_completed <- 0
    num_total <- job@arraysize
    while(TRUE) {
      info <- JobStatus(job)[,c("status","num")]
      info <- info[!stringr::str_detect(info$num, stringr::fixed("+")),,drop=FALSE] #ignore jobs like 1+

      #Merge new info with last info. Only keep the latest entries
      info <- rbind(info, last_info)
      info <- info[!duplicated(info$num),] #remove older entries
      last_info <- info
      
      #For each expected job in the array, figure out it's status
      current_state <- rep("PENDING",job@arraysize) ## keep outside??
      num_running <- 0
      num_failed <- 0
      num_outofmem <- 0
      num_cancelled <- 0
      num_timeout <- 0
      num_terminal <- 0
      
      if(nrow(info)==0) {
        #If we are extremely unlucky, this can suddenly happen overnight while waiting for jobs to be done
        cur_summary <- "Waiting for SLURM job to start"
      } else {
        
        #Parse each line of info
        for(i in 1:nrow(info)) {
          this_num <- info$num[i]
          if(this_num=="[0]") {
            #This is the entry for the batch array as a whole
            if(info$status[i] %in% slurm_failure_states()) {
              current_state[1:job@arraysize] <- info$status[i]
            }
          } else {
            
            #check if string is numeric
            if(grepl("^[0123456789]+$",this_num)) {
              
              #This is an individual index into the array
              slurm_index <- as.integer(this_num)
              if(slurm_index >= 0 && slurm_index < job@arraysize) {
                current_state[slurm_index + 1] <- info$status[i]
              }
              
            } else {
              #this might be a weird entry such as 35556921_0.+  ; not sure what to make of these
            }
          }
        }        
        
        #Count number of jobs of each type
        num_running <- floor(sum(current_state=="RUNNING"))
        num_completed <- floor(sum(current_state=="COMPLETED"))
        num_failed <- floor(sum(current_state=="FAILED"))
        num_outofmem <- floor(sum(current_state=="OUT_OF_MEM"))
        num_cancelled <- floor(sum(current_state=="CANCELLED"))
        num_timeout <- floor(sum(current_state=="TIMEOUT"))
        num_terminal <- floor(sum(current_state %in% slurm_terminal_states()))

        cur_summary <- paste0(
          job@pid," ",
          job@jobname,"   ",
          "To run: ",num_total,"   ",
          "Completed: ",num_completed,"   ",
          "Running: ",num_running,"   ",
          "Failed: ", num_failed,"   ",
          "Cancelled: ", num_cancelled,"   ",
          "Out-of-mem: ", num_outofmem,"   ",
          "Timeout: ", num_timeout
        )

      }
      
      cli::cli_progress_update(set = min(num_terminal, num_total))
      if(num_completed==num_total) {
        cli::cli_progress_done()
        cat("\nDone\n")
        break
      } else if(num_terminal==num_total) {
        cli::cli_progress_done()
        stop(
          paste0(
            "SLURM job ", job@pid, " (", job@jobname, ") finished with failed tasks: ",
            "Failed: ", num_failed, ", ",
            "Cancelled: ", num_cancelled, ", ",
            "Out-of-mem: ", num_outofmem, ", ",
            "Timeout: ", num_timeout
          ),
          call. = FALSE
        )
      } else {
        Sys.sleep(1)
      }
    }
  }
)

# job <- createSlurmJobFromExisting(35556921,2)
# JobStatus(job)
# WaitForJob(job)




#' @describeIn CancelJob SLURM job method.
#' @export
setMethod(
  f = "CancelJob",
  signature ="SlurmJob",
  definition = function(job) {
    #First check if it is running at all
    
    system(paste("scancel", job@pid))
  }
)





#' @describeIn JobLog SLURM job method.
#' @export
setMethod(
  f = "JobLog",
  signature ="SlurmJob",
  definition = function(job) {
    invisible()
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
  
  
  
  j <- createSlurmJobFromExisting("35772765", arraysize = 1)
  WaitForJob(j)
}
