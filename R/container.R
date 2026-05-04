################################################################################
################ Settings for the underlying Bascet installation ###############
################################################################################

###############################################
#' @export
setClass("BascetInstance", slots=list(
  bin="character",
  tempdir="character",
  prependCmd="character",
  containerMem="character",
  logLevel="character",
  logToFile="logical"
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
#' @param containerMem Amount of memory used by the container itself
#' @param logLevel ...
#' 
#' @return A Bascet instance
#' @export
BascetInstance <- function(
    bin, 
    tempdir, 
    prependCmd="",
    containerMem="0B",
    logLevel=c("info", "debug", "warn"),
    logToFile=TRUE
){
  #check arguments
  logLevel <- match.arg(logLevel)

  if(!is.null(tempdir) & !file.exists(tempdir)){
    stop(sprintf("temp directory %s does not exist", tempdir))
  }
  #cannot check other arguments; need to trust user

  stopifnot(is.valid.memsize(containerMem))
  
  #Generate instance
  new(
    "BascetInstance",
    bin=bin,
    tempdir=tempdir,
    prependCmd=prependCmd,
    containerMem=containerMem,
    logLevel=logLevel,
    logToFile=logToFile
  )
}


###############################################
#' Check that parameter is a valid bascet instance
#' @param x An object to test for BascetInstance class
#' @noRd
is.bascet.instance <- function(x) {
  stringr::str_detect(as.character(class(x)),"BascetInstance")
}

###############################################
# The default Bascet installation settings
# 
#bascetInstance.default <- BascetInstance(
#  bin="/home/mahogny/jupyter/bascet/target/debug/bascet",
#  tempdir="./"
#)



###############################################
#' Get default Bascet instance from global variable (bascetInstance.default)
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
#' 
#' @return A path to a temp directory that can be created. Must be removed when done
#' @export
GetBascetTempDir <- function(
    bascetInstance
){
  #Check arguments
  stopifnot(is.bascet.instance(bascetInstance))
  
  #Generate a tempfile if needed
  if(is.null(bascetInstance@tempdir)){
    tempfile()
  } else {
    tempfile(tmpdir=file.path(bascetInstance@tempdir))
  }
}







###############################################
#' Get a Bascet executable
#' It will be cached in the provided directory to avoid downloading it each the time the function is called
#' 
#' @param storeAt Directory to store the container in. Default is current directory but it is likely better to provide a single systems level directory #################### TODO
#' @param tempdir Default is to create a directory for temporary files in the current directory. Place it on a fast disk if possible
#' @param logLevel Log level for the Bascet instance (e.g. "info", "debug", "warn")
#'
#' @return A Bascet instance
#' @export
getBascetExecutable <- function(
    storeAt=getwd(),
    tempdir=NULL,
    logLevel="info"
){
  #Check arguments
  stopifnot(dir.exists(storeAt))

  if(is.null(tempdir)){
    tempdir <- "./temp"
    dir.create(tempdir, showWarnings = FALSE)
  } else {
    stopifnot(dir.exists(tempdir))
  }

  #Download file if needed
  #file_bascet_sif <- file.path(storeAt, "bascet.sif")
  #if(!file.exists(file_bascet_sif)) {
  #  print("No singularity image present; downloading")
  #  safeDownloadMD5("http://beagle.henlab.org/public/bascet/bascet.sif",file_bascet_sif)
  #} else {
  #  print(paste("Found existing Bascet singularity image:", file_bascet_sif))
  #}

  prependCmd <- ""  
  #paste("singularity run", file_bascet_sif," ")


  BascetInstance(
    bin="bascet",   ##### what about .exe ??? ################################  ###################### path?? ##################################
    tempdir=tempdir,
    prependCmd=prependCmd,
    containerMem="10GB",
    logLevel=logLevel
  )
}



###############################################
#' Get a Bascet binary for the current platform
#' It will be cached in the provided directory to avoid downloading it each the time the function is called
#'
#' @param storeAt Directory to store the binary in. Default is current directory but it is likely better to provide a single systems level directory
#' @param tempdir Default is to create a directory for temporary files in the current directory. Place it on a fast disk if possible
#' @param logLevel Log level for the Bascet instance (e.g. "info", "debug", "warn")
#' @param forceInstall Force download of the Bascet binary even if a cached binary exists
#'
#' @return A Bascet instance
#' @export
getBascetBinary <- function(
    storeAt=getwd(),
    tempdir=NULL,
    logLevel="info",
    forceInstall=FALSE
){
  #Check arguments
  stopifnot(dir.exists(storeAt))
  stopifnot(is.logical(forceInstall), length(forceInstall) == 1, !is.na(forceInstall))

  if(is.null(tempdir)){
    tempdir <- "./temp"
    dir.create(tempdir, showWarnings = FALSE)
  } else {
    stopifnot(dir.exists(tempdir))
  }

  root_url <- "http://beagle.henlab.org/public/bascet/bins"
  sysname <- tolower(Sys.info()[["sysname"]])
  file_local_bascet <- file.path(storeAt, "bascet")
  bin_path <- function(path) {
    normalizePath(path, mustWork=TRUE)
  }

  if(!forceInstall && file.exists(file_local_bascet)) {
    print(paste("Found existing Bascet binary:", file_local_bascet))
    if(sysname != "windows") {
      Sys.chmod(file_local_bascet, mode="0755")
    }
    return(BascetInstance(
      bin=bin_path(file_local_bascet),
      tempdir=tempdir,
      prependCmd="",
      containerMem="10GB",
      logLevel=logLevel
    ))
  }

  if(sysname == "linux") {
    bascet_bin <- "bascet-linux-x86_64"
  } else if(sysname == "darwin") {
    bascet_bin <- "bascet-macos-universal"
  } else if(sysname == "windows") {
    bascet_bin <- "bascet-windows-x86_64.exe"
  } else {
    stop(sprintf("Unsupported operating system for Bascet binary: %s", sysname))
  }

  file_bascet_bin <- file.path(storeAt, bascet_bin)
  url_bascet_bin <- paste(root_url, bascet_bin, sep="/")

  if(forceInstall || !file.exists(file_bascet_bin)) {
    if(forceInstall) {
      print("Force installing Bascet binary; downloading")
    } else {
      print("No Bascet binary present; downloading")
    }
    safeDownloadMD5(url_bascet_bin, file_bascet_bin)
  } else {
    print(paste("Found existing Bascet binary:", file_bascet_bin))
  }

  if(sysname != "windows") {
    Sys.chmod(file_bascet_bin, mode="0755")
  }

  BascetInstance(
    bin=bin_path(file_bascet_bin),
    tempdir=tempdir,
    prependCmd="",
    containerMem="10GB",
    logLevel=logLevel
  )
}





###############################################
#' Check if a Bascet instance works
#' 
#' @param bascetInstance Bascet instance
#' 
#' @return "ok" if the instance works; panic otherwise
#' @export
TestBascetInstance <- function(
    bascetInstance
) {
  #check arguments
  stopifnot(is.bascet.instance(bascetInstance))

  #prependCmd is a pre-built shell string (contains $PWD/%cd% expansion and
  #container map flags), so it must be interpreted by a shell.
  cmd <- paste(
    bascetInstance@prependCmd,
    shQuote(bascetInstance@bin),
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
  



###############################################
#' Download a file, check MD5 to ensure success. This assumes a file.md5 is stored on the server
#' 
#' @param url URL to the file to download
#' @param file Name of the file to download content to
#' 
#' @return Nothing; panics if the download fails
#' @noRd
safeDownloadMD5 <- function(
    url, 
    file
){
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
    stop("MD5 does not match the downloaded file")
  }
}



###############################################
#' Prepare Bascet command given arguments
#' 
#' @param bascetInstance bascetInstance
#' @param params List of parameters
#' 
#' @return Bascet command to run as a string
#' @noRd
assembleBascetCommand <- function(bascetInstance, params) {
  tor <- stringr::str_flatten(collapse = " ",
    c(
      bascetInstance@prependCmd,
      shQuote(bascetInstance@bin),
      if(bascetInstance@logToFile) "--log-mode=$BASCET_LOGFILE",
      paste0("--log-level=",bascetInstance@logLevel),
      params
    )
  )
#  print("///////////////")
#  print(tor)
#  print("///////////////")
  tor
}
