
################################################################################
################ Generics for running jobs #####################################
################################################################################



###############################################
#' Run a job
#'
#' Submit a job script with a runner backend.
#'
#' @param runner Runner object used to execute the job.
#' @param jobname Name for the submitted job.
#' @param bascetInstance Bascet instance that provides runtime settings.
#' @param cmd JobScript object to execute.
#' @param arraysize Number of array tasks.
#'
#' @return A job object, or a no-op job for synchronous/no-runner backends.
#' @export
setGeneric(
  name = "RunJob",
  def = function(runner, jobname, bascetInstance, cmd, arraysize) standardGeneric("RunJob")
)

###############################################
#' Wait for a job to finish
#'
#' @param job A job object.
#'
#' @return Invisibly returns when the job has completed.
#' @export
setGeneric(
  name = "WaitForJob",
  def = function(job) standardGeneric("WaitForJob")
)

###############################################
#' Cancel a job
#'
#' @param job A job object.
#'
#' @return Backend-specific cancellation result, usually invisible.
#' @export
setGeneric(
  name = "CancelJob",
  def = function(job) standardGeneric("CancelJob")
)

###############################################
#' Get job status
#'
#' @param job A job object.
#'
#' @return A data frame or backend-specific status object.
#' @export
setGeneric(
  name = "JobStatus",
  def = function(job) standardGeneric("JobStatus")
)

###############################################
#' Get job log
#'
#' @param job A job object.
#'
#' @return Backend-specific log content, usually invisible if no log is available.
#' @export
setGeneric(
  name = "JobLog",
  def = function(job) standardGeneric("JobLog")
)




################################################################################
########### No runner - for testing ############################################
################################################################################


###############################################
#' No-op runner
#'
#' Runner backend used for dry runs and debugging.
#'
#' @slot showScript Logical flag indicating whether scripts should be shown.
#' @export
setClass("NoRunner", slots=list(
  showScript="logical"
  )
) 



###############################################
#' @describeIn RunJob No-op runner method.
#' @export
setMethod(
  f = "RunJob",
  signature ="NoRunner",
  definition = function(runner, jobname, bascetInstance, cmd, arraysize) {  
    if(!is.jobscript(cmd)) {
      stop("NoRunner requires cmd to be a JobScript")
    }
    print("This is an attempt to start a job (but will not run)")
    new_no_job()
  }
)



###############################################
#' Create new no-runner instance, used for debugging
#' 
#' @param showScript For debugging: print script to run
#' 
#' @return A NoRunner instance
#' @export
NoRunner <- function(
    showScript=TRUE
){
  new("NoRunner", showScript=showScript)
}



################################################################################
########### No job - placeholder to return if a job was not needed #############
################################################################################

###############################################
#' No-op job
#'
#' Placeholder job used when no asynchronous process exists.
#'
#' @param object A NoJob object.
#' @slot evil Placeholder character slot to keep the S4 class concrete.
#' @export
setClass("NoJob", slots=list(
  evil="character"
)
) 

###############################################
#' Create an empty job. It is considered to have terminated from the start
#' 
#' @return A no-job
#' @noRd
new_no_job <- function() {
  new(
    "NoJob",
    evil="666" #might need a member to avoid virtual class errors
  )
}


###############################################
#Has possibility of ctrl+c; just keeps polling, possibly with a status indicator from log. or keep plotting log file
#' @describeIn WaitForJob No-op job method.
#' @export
setMethod(
  f = "WaitForJob",
  signature ="NoJob",
  definition = function(job) {
    invisible()
  }
)

###############################################
#' @describeIn CancelJob No-op job method.
#' @export
setMethod(
  f = "CancelJob",
  signature ="NoJob",
  definition = function(job) {
    invisible()
  }
)

###############################################
#' @describeIn JobStatus No-op job method.
#' @export
setMethod(
  f = "JobStatus",
  signature ="NoJob",
  definition = function(job) {
    invisible()
  }
)



###############################################
#' @describeIn NoJob Print a compact representation of a no-op job.
#' @export
setMethod(
  f = "show", 
  signature = "NoJob", 
  definition = function(object) {
    #Hide the content
    ""
  })



###############################################
#' Check that parameter is a valid runner
#' @param x An object to test for Runner class
#' @noRd
is.runner <- function(x) {
  stringr::str_detect(as.character(class(x)),"Runner")
}

################################################################################
########### Default runner #####################################################
################################################################################


###############################################
#' Get the current default runner
#' 
#' @return The current default runner
#' @export
GetDefaultBascetRunner <- function() {
  
  #Can create and return this if variable does not exist
  #bascetRunner.default <- LocalRunner(direct = TRUE, showScript=FALSE)
  
  
  bascetRunner.default
}

