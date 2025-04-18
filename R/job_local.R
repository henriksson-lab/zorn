

################################################################################
################ Things for running locally ####################################
################################################################################


#library(processx) ##### 


###############################################
#' @export
setClass("LocalRunner", slots=list(
  maxcpu="character",
  direct="logical",
  show_script="logical"
)
) 


###############################################
#' @export
setClass("LocalJob", slots=list(
  cmd="character",
  proc="ANY", #processx:process
  logfile="character",
  arraysize="numeric"
)
) 

###############################################
#' Create new local runner instance
#' 
#' @return TODO
#' @export
LocalRunner <- function(maxcpu="10", direct=FALSE, show_script=FALSE){
  new("LocalRunner", maxcpu=maxcpu, direct=direct, show_script=show_script)
}




###############################################
#' @export
setMethod(
  f = "RunJob",
  signature ="LocalRunner",
  definition = function(runner, jobname, cmd, arraysize) {  ######## todo likely remove witdata
    
    print("Starting local job")
    cmd <- stringr::str_flatten(c(
      cmd,
      ""
    ),"\n")
    
    print("----------- arg")
    if(runner@show_script){
      print("=============== final script start =================")
      cat(cmd)
      print("=============== final script end =================")
    }
    
    
    #### Concatenate all scripts over all TASK_ID
    all_cmd <- c()
    for(i in seq_len(arraysize)){
      this_cmd <- stringr::str_replace_all(cmd,stringr::fixed("$TASK_ID"),i-1)  #0... TASK_ID -1
    ##  onep <- processx::process$new(this_cmd)
      all_cmd <- c(all_cmd,this_cmd)
    }
    
    
    if(runner@show_script){
      print("=============== all script start =================")
      cat(all_cmd)
      print("=============== all script end =================")
    }
    
    
    tfile <- tempfile(pattern = jobname, fileext=".sh") #putting jobname here helps it show up in "ps"; but may cause issues if bad jobname given
    writeLines(con=tfile,c(
      "#!/bin/bash",
      stringr::str_flatten(all_cmd,"\n")
    ))
    
    if(runner@direct) {
      print("Running directly")  
      system(paste("bash", tfile))
      #print(tfile)
      file.remove(tfile) #assumes process has started. can we do better?
      
      #Return no job, as it is done already
      new_no_job()
    } else {
      print("Using separate process")
      write(tfile,paste("\n", "rm", tfile), append=TRUE)  ##this will make the script delete itself at the end ; or put command first?
      job <- new(
        "LocalJob",
        cmd=cmd,
        proc=processx::process$new("bash",tfile),
        logfile="foo",
        arraysize=arraysize
      )
      
      #Sys.sleep(1)
      #file.remove(tfile) #assumes process has started. can we do better? yes -- can put delete command into the tfile itself!! TODO
      
      #Return the job with PIDs set  
      job
    }
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

bascet_runner.default <- LocalRunner(direct = TRUE, show_script=FALSE)



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

