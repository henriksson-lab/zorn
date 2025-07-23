


# This is a way of expanding files from shell scripts, just as they are needed
# note that environmental variables can be at most 1M on osx
# 
# TASK_ID=2
# 
# TMPFILE=$(mktemp --directory)
# 
# if test $TASK_ID -eq 1; then
# echo -e "now it is 1" > ${TMPFILE}/cells.0.txt
# echo -e "now it is b" >> ${TMPFILE}/cells.0.txt
# fi
# 
# if test $TASK_ID -eq 2; then
# echo -e "now it is 2" > ${TMPFILE}/cells.0.txt
# fi
# 
# cat ${TMPFILE}/cells.0.txt


###############################################
#' Generate a shell script command to produce a file of list of strings
#' 
#' OSX limit for command line is 1048576 characters. limit to about 200k to be on the safe side
#' 
#' Possible to check limit by: getconf ARG_MAX
#' 
#' Note that the list need to start from [[1]]
#' 
#' This system only works for smaller lists; slurm complains if too large .sh-files are submitted
#' 
# shellscriptMakeFilesExpander <- function(tmpname, list_content, compare_to_var="TASK_ID") {
#   script_lines <- list()
#   script_lines[[1]] <- paste0(tmpname,"=$(mktemp)")
#   j <- 2
#   
#   rowlen <- 20000000 #temporarily disabled. this approach to making files is dangerous 
#   
#   for(i in 1:length(list_content)){
#     
#     this_str <- stringr::str_flatten(list_content[[i]], collapse = "\\n")
#     #this_str <- "123456789abcdefghijklmnopqrst123456789abcdefghijklmnopqrst123456789abcdefghijklmnopqrst123456789abcdefghijklmnopqrst123456789abcdefghijklmnopqrst"
#     
#     #Find borders for dividing the string
#     all_starts <- seq(1, stringr::str_length(this_str), by=rowlen)
#     itvs <- c(all_starts, stringr::str_length(this_str)+1)
#     all_ends <- itvs[-1] - 1
# 
#     
#     ##### bug!!! might split on a \??
#     
#     #Divide the string accordingly    
#     divided_str <- stringr::str_sub(this_str, start = all_starts, end = all_ends)
# 
#     #Produce an if-statement to generate a file based on TASK_ID or similar
#     script_lines[[j]] <- paste0("if test $",compare_to_var," -eq ",i-1,"; then")   #note, bash is indexing from 0
#     j <- j+1
# 
#     #For this list item, generate the echo statements for each text fragment
#     for(k in 1:length(divided_str)){
#       if(k==1) {
#         script_lines[[j]] <-  paste0(
#           "echo -e ",
#           "\"",divided_str[k], "\"",
#           " > ${",tmpname,"}")
#       } else {
#         script_lines[[j]] <-  paste0(
#           "echo -e ",
#           "\"",divided_str[k], "\"",
#           " >> ${",tmpname,"}")
#       }
#       j <- j+1
#     }
# 
#     script_lines[[j]] <- "fi"
#     j <- j+1
#     
#   }
#   
#   script_lines[[j]] <- "\n"
# 
#   tot <- stringr::str_flatten(script_lines, collapse = "\n")
#   tot
# }







###############################################
#' Generate a shell script command to produce a file of list of strings
#' 
#' The name of each file will be in VARIABLE[TASK_ID] where TASK_ID starts from 0
#' 
#' Note that there must be one file per task, or each task deletes its own file. This will result in race conditions
#' 
#' 
shellscriptMakeFilesExpander <- function(for_variable, list_content) {
  
  #Figure out the name of a tempfile to use.
  #Note that we cannot use the regular tempfile() call as these files would be deleted
  #when R is closed, possibly before SLURM or other tool has consumed them
  ts <- as.character(format(
    as.numeric(Sys.time())*1000, 
    scientific=FALSE
  ))
  #print(ts)
  
  #General directory to store these temporary files in 
  path_tmp <- file.path(".jobdata")
  dir.create(path_tmp, showWarnings = FALSE)
  
  #Create each temp file
  all_files <- c()
  for(i in 1:length(list_content)){
    one_tmp <- file.path(path_tmp, paste0(ts,"-",i))
    if(file.exists(one_tmp)) {
      stop(paste0("Attempting to create file ",one_tmp,", but it already exists"))
    }
    writeLines(con = one_tmp, text = list_content[[i]])
    all_files <- c(all_files, one_tmp)
  }
  
  #Create an array
  cmd1 <- shellscriptMakeBashArray(for_variable, all_files)
  
  #Delete the file upon exit
  cmd2 <- paste0(
    "trap \"rm -rf ${",for_variable,"[$TASK_ID]}\" EXIT"
  )

  #Return commands for shell script
  paste(
    cmd1,
    cmd2,
    sep = "\n"
  )
}











###############################################
#' from [a,b] to "a,b"
shellscriptMakeCommalist <- function(f) {
  stringr::str_flatten(f, collapse = ",")
}


###
#' Helper function: Create array of values in bash scripts
# example: myArray=("cat" "dog" "mouse" "frog")
shellscriptMakeBashArray <- function(variable, vals){
  paste0(
    variable,
    "=(",
    stringr::str_flatten(sprintf("\"%s\"",vals), collapse = " "),
    ")")  
}

if(FALSE){
  cat(shellscriptMakeBashArray("files_r2",c("a","b")))
}



###############################################
#' 
#' 
shellscriptSplitArrayIntoListRandomly <- function(arr, num_divide){
  set.seed(666)
  shard_assignment_for_cell <- sample(1:num_divide, length(arr), replace = TRUE)
  outlist <- list()
  for(i in 1:num_divide){
    outlist[[i]] <- arr[shard_assignment_for_cell==i]
  }
  outlist
}





###############################################
#' 
#' 
#shellscript_set_tempdir <- function(bascetInstance){
#  paste0("BASCET_TEMPDIR=",GetBascetTempDir(bascetInstance))
#}






###############################################
#' 
#' 
shellscriptMakeOneFileExpander <- function(tmpname, list_lines) {
  script_lines <- c(
    paste0(tmpname,"=$(mktemp)"),
    paste0(
      "echo -e ",
      "\"",stringr::str_flatten(list_lines, collapse = "\\n"), "\"",  
      " > ${",tmpname,"}")
  )
  tot <- stringr::str_flatten(script_lines, "\n")
  tot
}






###############################################
#' Create a piece of script to exit a job early if file exists
#' 
#' @return script 
shellscriptCancelJobIfFileExists <- function(fvar) {
  c(
    paste0(
      "if [ -f ",fvar," ]; then"),
    "  echo \"Skipping job as the output exists already\"",
    "  exit 0",
    "fi"
  )
}



if(FALSE){
  
  list_content <- list()
  list_content[[1]] <- c("1","3","3")
  list_content[[2]] <- c("5","3")
  
  tmpname <- "TMPFILE"
  
  cat(helper_makescript_files_expander(list_content, "TMPFILE"))
  
  cat(tot)
  
  
}
