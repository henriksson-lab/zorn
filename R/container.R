 


################################################################################
################ Settings for the underlying Bascet installation ###############
################################################################################

###############################################
#' @export
setClass("BascetInstance", slots=list(
  bin="character",
  tempdir="character",
  prependCmd="character"
)
) 

################################################################################
################ Functions #####################################################
################################################################################

###############################################
#' Create a new bascet instance.
#' For advanced users only
#' 
#' @param bin Name of the binary
#' @param tempdir Directory where to store temporary files
#' @param prependCmd Something to prepend to the command, to e.g. support container systems
#' @return A Bascet instance
#' @export
BascetInstance <- function(
    bin, 
    tempdir, 
    prependCmd=""
){
  if(!is.null(tempdir) & !file.exists(tempdir)){
    stop(sprintf("temp directory %s does not exist", tempdir))
  }
  new(
    "BascetInstance",
    bin=bin,
    tempdir=tempdir,
    prependCmd=prependCmd
  )
}



###############################################
# The default Bascet installation settings
# 
#bascetInstance.default <- BascetInstance(
#  bin="/home/mahogny/jupyter/bascet/target/debug/bascet",
#  tempdir="./"
#)



###############################################
#' Get Bascet instance from global variable
#' 
#' @return A Bascet instance
#' @export
GetDefaultBascetInstance <- function(){
  bascetInstance.default
}



###############################################
#' Get a temp directory to use; need to be created
#' 
#' @param bascetInstance A Bascet instance
#' @return A path to a temp directory that can be created. Must be removed when done
#' @export
GetBascetTempDir <- function(
    bascetInstance
){
  if(is.null(bascetInstance@tempdir)){
    tempfile()
  } else {
    tempfile(tmpdir=file.path(bascetInstance@tempdir))
  }
}





###############################################
#' Get a Bascet image (singularity or docker). 
#' It will be cached in the provided directory to avoid downloading it all the time
#' 
#' @return A Bascet instance
#' @export
getBascetSingularityImage <- function(
    storeAt=getwd(), 
    tempdir=NULL
){
  
  file_bascet_sif <- file.path(storeAt, "bascet.sif")
  
  if(!file.exists(file_bascet_sif)) {
    print("No singularity image present; downloading")
 
    safeDownloadMD5("http://beagle.henlab.org/public/bascet/bascet.sif",file_bascet_sif)

  } else {
    print(paste("Found existing Bascet singularity image:", file_bascet_sif))
  }
  
  prependCmd <- paste("singularity run", file_bascet_sif," ")
  
  if(is.null(tempdir)){
    tempdir <- "./"
  }

  BascetInstance(
    bin="bascet",
    tempdir=tempdir,
    prependCmd=prependCmd
  ) 
}

if(FALSE){
  getBascetDockerImage(storeAt = "~/mystore/temp")
  getBascetSingularityImage(storeAt = "~/mystore/temp")
}


###############################################
#' Get and install a Bascet docker image.
#' 
#' @return A Bascet instance
#' @export
getBascetDockerImage <- function(
    storeAt=getwd(),
    tempdir=NULL,
    forceInstall=FALSE,
    verbose=FALSE
){
  
  # docker image inspect busybox:latest >/dev/null 2>&1 && echo yes || echo no
  
  #check if docker is installed
  ret <- system("docker ps", ignore.stdout = !verbose, ignore.stderr = !verbose)
  if(ret==0) {
    ret <- system("docker image inspect henriksson-lab/bascet:latest", ignore.stdout = !verbose, ignore.stderr = !verbose)
    
    if(ret!=0 || forceInstall) {
      print("No docker image present; downloading")
      
      file_bascet_image <- file.path(storeAt, "bascet.tar")
      safeDownloadMD5("http://beagle.henlab.org/public/bascet/bascet.tar",file_bascet_image)

      print("Loading image into Docker")
      system(paste("docker load -i ", file_bascet_image))
      
      print(paste("The large image at",file_bascet_image,"can now be removed if the installation worked. You can otherwise try to install it manually using Docker"))
    } else {
      print(paste("Found existing Bascet Docker image"))
    }
    
    prependCmd <- paste("docker run henriksson-lab/bascet ")
    
    if(is.null(tempdir)){
      tempdir <- "./"
    }
    
    BascetInstance(
      bin="bascet",
      tempdir=tempdir,
      prependCmd=prependCmd
    ) 
    
  } else {
    stop("Docker is not installed or cannot be run")
  } 
}



###############################################
#' Remove current Bascet docker image 
#' 
#' @export
removeBascetDockerImage <- function(){
  system("docker image rm -f henriksson-lab/bascet")
}



###############################################
#' Check if a Bascet instance works
#' 
#' @param bascetInstance Bascet instance
#' @return "ok" if the instance works; panic otherwise
#' @export
TestBascetInstance <- function(
    bascetInstance
) {
  
  cmd <- paste(
    bascetInstance@prependCmd,
    bascetInstance@bin,
    "-V"
  )

  ret <- system(cmd, intern = TRUE)

  #Print version number (if it works)
  print(ret)
  
  if(stringr::str_detect(ret[1], "bascet")){
    "ok"
  } else {
    stop("Could not invoke Bascet")
  }
}
  



if(FALSE){
  bascet_inst <- getBascetSingularityImage("/data/henlab/temp/")
  TestBascetInstance(bascet_inst)
  
  #singularity run bascet_0.01.sif skesa
  #singularity exec lolcow_latest.sif cowsay moo
}



###############################################
#' Download a file, check MD5 to ensure success
#' 
#' @return nothing; panics if the download fails
safeDownloadMD5 <- function(url, file){
  url_md5 <- paste0(url,".md5")
  file_md5 <- paste0(file,".md5")
  
  f <- RCurl::CFILE(file_md5, mode="wb")
  a <- RCurl::curlPerform(url = url_md5, writedata = f@ref, noprogress=FALSE)
  RCurl::close(f)
  
  line_md5 <- readLines(file_md5)
  
  if(file.exists(file_md5)){
    file.remove(file_md5)
  }
  
  if(length(line_md5)>1) {
    print(line_md5)
    stop("Error in reading the MD5 file")
  }
  
  prev_md5 <- stringr::str_split_fixed(line_md5, " ",2)[1]
  print(paste("Got previous MD5 value to compare against:", prev_md5))
  #return(a)
  
  
  
  f <- RCurl::CFILE(file, mode="wb")
  a <- RCurl::curlPerform(url = url, writedata = f@ref, noprogress=FALSE)
  RCurl::close(f)
  
  print("Computing MD5 for downloaded file")
  new_md5 <- unname(tools::md5sum(file))
  print(paste("MD5 is:", new_md5))
  
  if(new_md5==prev_md5) {
    print("MD5 matches")
  } else {
    if(file.exists(file)){
      file.remove(file)
    }
    stop("MD5 does not match")
  }
  #return(TRUE)
}


