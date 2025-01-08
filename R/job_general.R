
################################################################################
################ Generics for running jobs #####################################
################################################################################



setGeneric(
  name = "RunJob",
  def = function(runner, jobname, cmd, arraysize) standardGeneric("RunJob")
)

setGeneric(
  name = "WaitForJob",
  def = function(job) standardGeneric("WaitForJob")
)

setGeneric(
  name = "CancelJob",
  def = function(job) standardGeneric("CancelJob")
)

setGeneric(
  name = "JobStatus",
  def = function(job) standardGeneric("JobStatus")
)

setGeneric(
  name = "JobLog",
  def = function(job) standardGeneric("JobLog")
)


