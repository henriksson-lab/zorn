
################################################################################
################ Generics for running jobs #####################################
################################################################################



#' @export
setGeneric(
  name = "RunJob",
  def = function(runner, jobname, cmd, arraysize) standardGeneric("RunJob")
)

#' @export
setGeneric(
  name = "WaitForJob",
  def = function(job) standardGeneric("WaitForJob")
)

#' @export
setGeneric(
  name = "CancelJob",
  def = function(job) standardGeneric("CancelJob")
)

#' @export
setGeneric(
  name = "JobStatus",
  def = function(job) standardGeneric("JobStatus")
)

#' @export
setGeneric(
  name = "JobLog",
  def = function(job) standardGeneric("JobLog")
)


