


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
#' 
#' 
shellscript_make_files_expander <- function(tmpname, list_content, compare_to_var="TASK_ID") {
  script_lines <- c(
    paste0(tmpname,"=$(mktemp)")
  )
  for(i in 1:length(list_content)){
    script_lines <- c(
      script_lines,
      paste0("if test $",compare_to_var," -eq ",i-1,"; then"),   #note, indexing from 0
      paste0(
        "echo -e ",
        "\"",stringr::str_flatten(list_content[[i]], collapse = "\\n"), "\"",  ##list has to start from 1 !!
        " > ${",tmpname,"}"),
      "fi"
    )
  }
  tot <- stringr::str_flatten(script_lines, "\n")
  tot
}

###############################################
#' from [a,b] to "a,b"
shellscript_make_commalist <- function(f) {
  stringr::str_flatten(f, collapse = ",")
}


###
#' Helper function: Create array of values in bash scripts
# example: myArray=("cat" "dog" "mouse" "frog")
shellscript_make_bash_array <- function(variable, vals){
  paste0(
    variable,
    "=(",
    stringr::str_flatten(sprintf("\"%s\"",vals), collapse = " "),
    ")")  
}

#cat(make_bash_array("files_r2",c("a","b")))



###############################################
#' 
#' 
shellscript_split_arr_into_list_randomly <- function(arr, num_divide){
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
shellscript_set_tempdir <- function(bascet_instance){
  paste0("BASCET_TEMPDIR=",GetBascetTempDir(bascet_instance))
}






###############################################
#' 
#' 
shellscript_make_one_file_expander <- function(tmpname, list_lines) {
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
helper_cancel_job_if_file_exists <- function(fvar) {
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
