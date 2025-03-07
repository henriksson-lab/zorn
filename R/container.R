 


################################################################################
################ Settings for the underlying Bascet installation ###############
################################################################################

#' @export
setClass("BascetInstance", slots=list(
  bin="character",
  tempdir="character",
  prepend_cmd="character"
)
) 


###############################################
#' Create a new bascet instance
#' @return TODO
#' @export
BascetInstance <- function(bin, tempdir, prepend_cmd=""){
  if(!is.null(tempdir) & !file.exists(tempdir)){
    stop(sprintf("temp directory %s does not exist", tempdir))
  }
  new(
    "BascetInstance",
    bin=bin,
    tempdir=tempdir,
    prepend_cmd=prepend_cmd
  )
}



###############################################
#' The default Bascet installation settings
#' @return TODO
#' @export
bascet_instance.default <- BascetInstance(
  bin="/home/mahogny/jupyter/bascet/target/debug/robert",
  tempdir="/data/henlab/bascet_temp"
)

#bascet_instance.default <- BascetInstance(bin="bascet")




###############################################
#' Get a temp directory to use; need to be created
#' @return TODO
#' @export
GetBascetTempDir <- function(bascet_instance){
  if(is.null(bascet_instance@tempdir)){
    tempfile()
  } else {
    tempfile(tmpdir=file.path(bascet_instance@tempdir))
  }
}






###############################################
#' Get a Bascet image (singularity or docker). 
#' It will be cached in the provided directory to avoid downloading it all the time
#' 
#' TODO Should be managed by Bascet mapshard system, with automatic input conversion. unaligned file should be made temp and removed
#' 
#' @return TODO
#' @export
getBascetImageInstance <- function(store_at="./", tempdir=NULL) {
  
  
  file_bascet_sif <- file.path(store_at, "bascet_0.01.sif")

  if(!file.exists(file_bascet_sif)) {
    print("No singularity image present; downloading")
    ret <- system("singularity pull --arch amd64 library://lmc297/bascet/bascet:0.01")

    if(ret==127){
      #TODO: check for errors

      stop("Failed to download singularity")      
    }
        
  } else {
    print(paste("Found existing Bascet singularity image:", file_bascet_sif))
  }
 
  prepend_cmd <- paste("singularity run ", file_bascet_sif," ")
  
  if(is.null(tempdir)){
    tempdir <- tempdir()
#    tempdir <- "/data/henlab/bascet_temp" ##TODO better place?
  }
  
  BascetInstance(
    bin="bascet",
    tempdir=tempdir,
    prepend_cmd=prepend_cmd
  )
}

if(FALSE){
  inst <- getBascetImageInstance()
  
  inst

  #singularity run bascet_0.01.sif skesa
  #singularity exec lolcow_latest.sif cowsay moo

  #might exec be faster? is miniconda the issue?
    
}






