
################################################################################
################ Generics for running jobs #####################################
################################################################################



setGeneric(
  name = "RunJob",
  def = function(runner, jobname, withdata, cmd, arraysize) standardGeneric("RunJob")
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

### Figure out how many chards


setClass("BascetInstance", slots=list(
  bin="character"
  )
) 

BascetInstance <- function(bin){
  new(
    "BascetInstance",
    bin=bin
  )
}

#bascet_instance.default <- BascetInstance(bin="bascet")
bascet_instance.default <- BascetInstance(bin="/home/mahogny/jupyter/bascet/target/debug/robert")



