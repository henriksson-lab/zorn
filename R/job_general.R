
################################################################################
################ Generics for running jobs #####################################
################################################################################



###############################################
#' @export
setGeneric(
  name = "RunJob",
  def = function(runner, jobname, cmd, arraysize) standardGeneric("RunJob")
)

###############################################
#' @export
setGeneric(
  name = "WaitForJob",
  def = function(job) standardGeneric("WaitForJob")
)

###############################################
#' @export
setGeneric(
  name = "CancelJob",
  def = function(job) standardGeneric("CancelJob")
)

###############################################
#' @export
setGeneric(
  name = "JobStatus",
  def = function(job) standardGeneric("JobStatus")
)

###############################################
#' @export
setGeneric(
  name = "JobLog",
  def = function(job) standardGeneric("JobLog")
)




################################################################################
########### No job - placeholder to return if a job was not needed #############
################################################################################




#' @export
setClass("NoRunner", slots=list(
  show_script="logical"
  )
) 


###############################################
#' Create new No-runner instance. This system is for testing only
#' 
#' @return A no-runner
#' @export
NoRunner <- function(show_script=FALSE){
  new("NoRunner",show_script=show_script)
}




#' @export
setMethod(
  f = "RunJob",
  signature ="NoRunner",
  definition = function(runner, jobname, cmd, arraysize) {  
    print("Attempt to start a job (will not run)")
    new_no_job()
  }
)




#' @export
setClass("NoJob", slots=list(
  evil="character"
)
) 

###############################################
#' Create an empty job. It is considered to have terminated from the start
#' 
#' @return A no-job
new_no_job <- function() {
  new(
    "NoJob",
    evil="666" #might need a member to avoid virtual class errors
  )
}





#Has possibility of ctrl+c; just keeps polling, possibly with a status indicator from log. or keep plotting log file
#' @export
setMethod(
  f = "WaitForJob",
  signature ="NoJob",
  definition = function(job) {
  }
)





#' @export
setMethod(
  f = "CancelJob",
  signature ="NoJob",
  definition = function(job) {
  }
)




#' @export
setMethod(
  f = "JobStatus",
  signature ="NoJob",
  definition = function(job) {
    NULL
  }
)
