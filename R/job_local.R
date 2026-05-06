

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
  proc="ANY",
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
    direct=NULL, 
    showScript=NULL
#    logToFile=NULL
){
  #Bring in defaults
  if(!is.null(settings) & is.null(ncpu)) {
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
      print(paste("==== mem ",mem))
    }
  } else {
    #Check that memory can be parsed and is some bare minimum
    stopifnot(parse_size_string(mem)>1e8)
  }
  
  if(is.null(showScript)) {
    showScript <- FALSE
  }
  
  if(is.null(direct)) {
    direct <- TRUE
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

    if(!is.jobscript(cmd)) {
      stop("LocalRunner requires cmd to be a JobScript")
    }

    if(!runner@direct) {
      stop("LocalRunner(direct = FALSE) is no longer supported")
    }

    if(runner@showScript) {
      print("=============== local JobScript start =================")
      str(cmd)
      print("=============== local JobScript end =================")
    }

    print("Running directly")
    tmproot <- GetBascetTempDir(bascetInstance)
    for(i in seq_len(arraysize)) {
      task_id <- i - 1L
      path_tempdir <- file.path(tmproot, as.character(task_id))
      dir.create(path_tempdir, recursive = TRUE, showWarnings = FALSE)

      task_env <- list(BASCET_TEMPDIR = path_tempdir)
      if(bascetInstance@logToFile) {
        dir.create("logs", showWarnings = FALSE)
        task_env$BASCET_LOGFILE <- paste0("logs/", jobname, ".", task_id, ".log")
      }

      tryCatch(
        runJobScriptLocal(cmd, task_id = task_id, env = task_env),
        error = function(e) {
          stop(
            paste0(
              "Local array task ",
              task_id,
              " failed while running ",
              jobname,
              ": ",
              conditionMessage(e)
            ),
            call. = FALSE
          )
        }
      )
    }

    new_no_job()
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
