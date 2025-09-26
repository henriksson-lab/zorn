
################################################################################
################ Generics for running jobs #####################################
################################################################################



###############################################
#' @export
setGeneric(
  name = "RunJob",
  def = function(runner, jobname, bascetInstance, cmd, arraysize) standardGeneric("RunJob")
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
########### No runner - for testing ############################################
################################################################################


###############################################
#' @export
setClass("NoRunner", slots=list(
  showScript="logical"
  )
) 



###############################################
#' @export
setMethod(
  f = "RunJob",
  signature ="NoRunner",
  definition = function(runner, jobname, bascetInstance, cmd, arraysize) {  
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
NoRunner <- function(showScript=TRUE){
  new("NoRunner", showScript=showScript)
}



################################################################################
########### No job - placeholder to return if a job was not needed #############
################################################################################

###############################################
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


###############################################
#Has possibility of ctrl+c; just keeps polling, possibly with a status indicator from log. or keep plotting log file
#' @export
setMethod(
  f = "WaitForJob",
  signature ="NoJob",
  definition = function(job) {
    invisible()
  }
)

###############################################
#' @export
setMethod(
  f = "CancelJob",
  signature ="NoJob",
  definition = function(job) {
    invisible()
  }
)

###############################################
#' @export
setMethod(
  f = "JobStatus",
  signature ="NoJob",
  definition = function(job) {
    invisible()
  }
)



###############################################
#' @export
setMethod(
  f = "show", 
  signature = "NoJob", 
  definition = function(object) {
    #Hide the content
    ""
  })



################################################################################
########### Default runner #####################################################
################################################################################


###############################################
#' Get the current default runner
#' 
#' @return The current default runner
#' @export
GetDefaultBascetRunner <- function() {
  bascetRunner.default
}

