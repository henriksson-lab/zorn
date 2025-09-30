


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
#' Check that parameter is a valid environment variable name
is.valid.env.variable <- function(x) {
  #TODO check that the name has no weird symbols
  is.character(x)
}


###############################################
#' Generate a shell script command to produce a file of list of strings
#' 
#' The name of each file will be in VARIABLE[TASK_ID] where TASK_ID starts from 0
#' 
#' Note that there must be one file per task, or each task deletes its own file. This will result in race conditions
#' 
#' @param for_variable Environment variable with array of files having the content
#' @param list_content Content for each file
#' 
#' @return Part of a shell script
shellscriptMakeFilesExpander <- function(
    for_variable, 
    list_content
) {
  #check arguments
  stopifnot(is.valid.env.variable(for_variable))
  stopifnot(is.list(list_content))
  
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
#' Helper function, taking a list of elements such as [a,b], and returning "a,b"
#' 
#' @param f List of elements
#' 
#' @return Formatted list
shellscriptMakeCommalist <- function(
    f
) {
  stringr::str_flatten(f, collapse = ",")
}


###############################################
#' Helper function: Create array of values in bash scripts
#' example: myArray=("cat" "dog" "mouse" "frog")
#' 
#' @param variable Environment variable to store the list in
#' @param vals Array of strings to store into the environment variable
#' 
#' @return A string similar to: myArray=("cat" "dog" "mouse" "frog")
shellscriptMakeBashArray <- function(
    variable, 
    vals
){
  #check arguments
  stopifnot(is.valid.env.variable(variable))
  
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
#' Helper function to take an array of elements and split it randomly into a number
#' of subset lists
#' 
#' @param arr Array of elements to subdivide
#' @param num_divide Number of subdivions
#' 
#' @return List of arrays
shellscriptSplitArrayIntoListRandomly <- function(
    arr, 
    num_divide
){
  set.seed(666)
  shard_assignment_for_cell <- sample(1:num_divide, length(arr), replace = TRUE)
  outlist <- list()
  for(i in 1:num_divide){
    outlist[[i]] <- arr[shard_assignment_for_cell==i]
  }
  outlist
}







###############################################
#' Helper function that takes content of a list and generates a BASH script 
#' that stores the content in a temporary file during execution
#' 
#' @param tmpname Name of environment variable in which the name of the file will be stored
#' @param list_lines Content to write to the file
#' 
#' @return BASH script content for the expander
shellscriptMakeOneFileExpander <- function(
    tmpname, 
    list_lines
) {
  #check arguments
  stopifnot(is.valid.env.variable(tmpname))
  
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
#' @param fvar Name of file
#' 
#' @return BASH script content 
shellscriptCancelJobIfFileExists <- function(
    fvar
) {
  #check arguments. cannot yet check if the file exists
  stopifnot(is.character(fvar))
  
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