 


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
#' The default Bascet installation settings
#' 
#' @return A Bascet instance
#' @export
bascetInstance.default <- BascetInstance(
  bin="/home/mahogny/jupyter/bascet/target/debug/bascet",
  tempdir="./"
)



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
 
    if(download.file("http://beagle.henlab.org/public/bascet/bascet.sif",file_bascet_sif,cacheOK = FALSE, mode = "wb")!=0){
      stop("Failed to download singularity image")      
    }
    
#    ret <- system("singularity pull --arch amd64 library://lmc297/bascet/bascet:0.01")
#    if(ret==127){
#      #TODO: check for errors
#      stop("Failed to download singularity")      
#    }
    
  } else {
    print(paste("Found existing Bascet singularity image:", file_bascet_sif))
  }
  
  prependCmd <- paste("singularity run", file_bascet_sif," ")
  
  if(is.null(tempdir)){
    tempdir <- tempdir()
  }
  
  BascetInstance(
    bin="bascet",
    tempdir=tempdir,
    prependCmd=prependCmd
  ) 
}




###############################################
#' Get and install a Bascet docker image.
#' 
#' @return A Bascet instance
#' @export
getBascetDockerImage <- function(
    storeAt=getwd(),
    tempdir=NULL
){
  
  # docker image inspect busybox:latest >/dev/null 2>&1 && echo yes || echo no

  #check if docker is installed
  ret <- system("docker ps")
  if(ret==0) {
    
    
    ret <- system("docker image inspect henriksson-lab/bascet:latest")
    
    if(ret!=0) {
      print("No docker image present; downloading")
      
      
      file_bascet_image <- file.path(storeAt, "bascet.tar")
      
      options(timeout = 60*60*5) #timeout in seconds
      
      if(download.file("http://beagle.henlab.org/public/bascet/bascet.tar",file_bascet_image, cacheOK = FALSE, mode = "wb")!=0){
        stop("Failed to download docker image")      
      }
      
      system(paste("docker load -i ", file_bascet_image))
      
      print(paste("The huge image at",file_bascet_image,"can now be removed if the installation worked. You can otherwise try to install it manually using Docker"))
      
    } else {
      print(paste("Found existing Bascet Docker image"))
    }
    
    prependCmd <- paste("docker run henriksson-lab/bascet ")
    
    if(is.null(tempdir)){
      tempdir <- tempdir()
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


