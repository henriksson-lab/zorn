################################################################################
################ Settings for the underlying Bascet installation ###############
################################################################################

###############################################
#' Bascet instance settings
#'
#' Runtime settings for the underlying Bascet installation.
#'
#' @slot bin Name or path of the Bascet binary.
#' @slot tempdir Directory for temporary files.
#' @slot prependCmd Command prefix, for example for container execution.
#' @slot containerMem Memory reserved for container overhead.
#' @slot logLevel Log level passed to Bascet.
#' @slot logToFile Whether Bascet logs are written to files.
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
#' @param logToFile Whether Bascet logs are written to files
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
  if(!is.null(tempdir)) {
    tempdir <- normalizeExistingDir(tempdir)
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
#' Get a Bascet binary from the target from a locally built Bascet repository
#' 
#' @param devdir Path to a local Bascet source checkout
#' @param tempdir Default is to create a directory for temporary files in the current directory. Place it on a fast disk if possible
#' @param logLevel Log level for the Bascet instance (e.g. "info", "debug", "warn")
#' @param targetType What target type to load
#' @param containerMem Amount of memory used by the container itself
#'
#' @return A Bascet instance
#' @export
getBascetDevDir <- function(
    devdir,
    tempdir=NULL,
    logLevel="info",
    targetType="release",
    containerMem="2GB"
){
  #Check arguments
  stopifnot(dir.exists(devdir))
  stopifnot(is.valid.memsize(containerMem))

  if(is.null(tempdir)){
    tempdir <- normalizePath("./temp", mustWork = FALSE)
    dir.create(tempdir, showWarnings = FALSE)
  } else {
    stopifnot(dir.exists(tempdir))
  }

  #Pull out target binary. Should be a file
  bascet_exe <- normalizePath(file.path(devdir, "target", targetType, "bascet"))
  stopifnot(file.exists(bascet_exe) & !dir.exists(bascet_exe))

  BascetInstance(
    bin=bascet_exe,
    tempdir=tempdir,
    prependCmd="",
    containerMem=containerMem,
    logLevel=logLevel
  )
}



###############################################
#' Get the default directory for storing Bascet binaries
#'
#' @return A directory path
#' @export
defaultBascetBinDir <- function() {
  candidate_dirs <- if(.Platform$OS.type == "windows") {
    c(
      file.path(Sys.getenv("USERPROFILE"), "Downloads"),
      file.path(path.expand("~"), "Downloads"),
      Sys.getenv("USERPROFILE"),
      path.expand("~")
    )
  } else {
    c(
      file.path(path.expand("~"), "Downloads"),
      path.expand("~")
    )
  }

  candidate_dirs <- candidate_dirs[nzchar(candidate_dirs)]
  candidate_dirs <- path.expand(candidate_dirs)
  candidate_dirs <- candidate_dirs[dir.exists(candidate_dirs)]

  if(length(candidate_dirs) == 0) {
    stop("Could not find a default directory for storing Bascet binaries")
  }

  candidate_dirs[[1]]
}



###############################################
#' Get a Bascet binary for the current platform from the development server
#' It will be cached in the provided directory to avoid downloading it each the time the function is called
#'
#' @param storeAt Directory to store the binary in. If NULL, uses \code{\link{defaultBascetBinDir}}
#' @param tempdir Default is to create a directory for temporary files in the current directory. Place it on a fast disk if possible
#' @param logLevel Log level for the Bascet instance (e.g. "info", "debug", "warn")
#' @param forceInstall Force download of the Bascet binary even if a cached binary exists
#' @param containerMem Amount of memory used by the container itself
#'
#' @return A Bascet instance
#' @export
getBascetBinaryDev <- function(
    storeAt=NULL,
    tempdir=NULL,
    logLevel="info",
    forceInstall=FALSE,
    containerMem="2GB"
){
  #Check arguments
  if(is.null(storeAt)) {
    storeAt <- defaultBascetBinDir()
  }
  stopifnot(dir.exists(storeAt))
  stopifnot(is.logical(forceInstall), length(forceInstall) == 1, !is.na(forceInstall))
  stopifnot(is.valid.memsize(containerMem))

  if(is.null(tempdir)){
    tempdir <- "./temp"
    dir.create(tempdir, showWarnings = FALSE)
  } else {
    stopifnot(dir.exists(tempdir))
  }
  tempdir <- normalizeExistingDir(tempdir)

  root_url <- "http://beagle.henlab.org/public/bascet/bins"
  sysname <- tolower(Sys.info()[["sysname"]])
  file_local_bascet <- file.path(storeAt, "bascet")
  bin_path <- function(path) {
    normalizePath(path, mustWork=TRUE)
  }

  if(!forceInstall && file.exists(file_local_bascet) && !dir.exists(file_local_bascet)) {
    print(paste("Found existing Bascet binary:", file_local_bascet))
    if(sysname != "windows") {
      Sys.chmod(file_local_bascet, mode="0755")
    }
    return(BascetInstance(
      bin=bin_path(file_local_bascet),
      tempdir=tempdir,
      prependCmd="",
      containerMem=containerMem,
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
    containerMem=containerMem,
    logLevel=logLevel
  )
}


###############################################
#' Get a Bascet binary for the current platform from GitHub releases
#' It will be cached in the provided directory to avoid downloading it each the time the function is called
#'
#' @param storeAt Directory to store the binary in. If NULL, uses \code{\link{defaultBascetBinDir}}
#' @param tempdir Default is to create a directory for temporary files in the current directory. Place it on a fast disk if possible
#' @param logLevel Log level for the Bascet instance (e.g. "info", "debug", "warn")
#' @param forceInstall Force download of the Bascet binary even if a cached binary exists
#' @param containerMem Amount of memory used by the container itself
#' @param repo GitHub repository in OWNER/REPO form
#'
#' @return A Bascet instance
#' @export
getBascetBinary <- function(
    storeAt=NULL,
    tempdir=NULL,
    logLevel="info",
    forceInstall=FALSE,
    containerMem="2GB",
    repo="henriksson-lab/bascet"
){
  #Check arguments
  if(is.null(storeAt)) {
    storeAt <- defaultBascetBinDir()
  }
  stopifnot(dir.exists(storeAt))
  stopifnot(is.logical(forceInstall), length(forceInstall) == 1, !is.na(forceInstall))
  stopifnot(is.valid.memsize(containerMem))
  stopifnot(is.character(repo), length(repo) == 1, !is.na(repo), nzchar(repo))

  if(is.null(tempdir)){
    tempdir <- "./temp"
    dir.create(tempdir, showWarnings = FALSE)
  } else {
    stopifnot(dir.exists(tempdir))
  }
  tempdir <- normalizeExistingDir(tempdir)

  sysname <- tolower(Sys.info()[["sysname"]])
  machine <- tolower(Sys.info()[["machine"]])
  file_local_bascet <- file.path(storeAt, if(sysname == "windows") "bascet.exe" else "bascet")
  bin_path <- function(path) {
    normalizePath(path, mustWork=TRUE)
  }

  if(!forceInstall && file.exists(file_local_bascet) && !dir.exists(file_local_bascet)) {
    print(paste("Found existing Bascet binary:", file_local_bascet))
    if(sysname != "windows") {
      Sys.chmod(file_local_bascet, mode="0755")
    }
    return(BascetInstance(
      bin=bin_path(file_local_bascet),
      tempdir=tempdir,
      prependCmd="",
      containerMem=containerMem,
      logLevel=logLevel
    ))
  }

  if(sysname == "linux") {
    bascet_archive <- "bascet-linux-x86_64.tar.gz"
  } else if(sysname == "darwin" && machine %in% c("arm64", "aarch64")) {
    bascet_archive <- "bascet-macos-aarch64.tar.gz"
  } else if(sysname == "darwin") {
    bascet_archive <- "bascet-macos-x86_64.tar.gz"
  } else if(sysname == "windows") {
    bascet_archive <- "bascet-windows-x86_64.zip"
  } else {
    stop(sprintf("Unsupported operating system for Bascet binary: %s", sysname))
  }

  release_url <- sprintf("https://github.com/%s/releases/latest/download", repo)
  file_archive <- file.path(storeAt, bascet_archive)
  url_archive <- paste(release_url, bascet_archive, sep="/")
  url_checksums <- paste(release_url, "checksums.txt", sep="/")

  if(forceInstall || !file.exists(file_archive)) {
    if(forceInstall) {
      print("Force installing Bascet binary; downloading")
    } else {
      print("No Bascet binary present; downloading")
    }
    safeDownloadSHA256(url_archive, file_archive, url_checksums)
  } else {
    print(paste("Found existing Bascet archive:", file_archive))
  }

  extract_dir <- file.path(storeAt, tools::file_path_sans_ext(tools::file_path_sans_ext(bascet_archive)))
  if(dir.exists(extract_dir)) {
    unlink(extract_dir, recursive=TRUE)
  }

  if(grepl("\\.zip$", bascet_archive)) {
    utils::unzip(file_archive, exdir=storeAt)
  } else {
    utils::untar(file_archive, exdir=storeAt)
  }

  extracted_bascet <- file.path(extract_dir, if(sysname == "windows") "bascet.exe" else "bascet")
  stopifnot(file.exists(extracted_bascet) && !dir.exists(extracted_bascet))
  file.copy(extracted_bascet, file_local_bascet, overwrite=TRUE)

  if(sysname != "windows") {
    Sys.chmod(file_local_bascet, mode="0755")
  }

  BascetInstance(
    bin=bin_path(file_local_bascet),
    tempdir=tempdir,
    prependCmd="",
    containerMem=containerMem,
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

  cmd <- buildBascetInvocation(bascetInstance, "-V")
  ret <- runSystem2Local(cmd$command, cmd$args, capture = TRUE)

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
#' Download a file and check SHA-256 against a GitHub release checksums.txt file.
#'
#' @param url URL to the file to download
#' @param file Name of the file to download content to
#' @param checksumsUrl URL to checksums.txt
#'
#' @return Nothing; panics if the download fails
#' @noRd
safeDownloadSHA256 <- function(
    url,
    file,
    checksumsUrl
){
  file_checksums <- paste0(file, ".checksums.txt")

  f <- RCurl::CFILE(file_checksums, mode="wb")
  a <- RCurl::curlPerform(url = checksumsUrl, writedata = f@ref, noprogress=FALSE)
  RCurl::close(f)

  checksum_lines <- readLines(file_checksums)

  if(file.exists(file_checksums)){
    file.remove(file_checksums)
  }

  filename <- basename(file)
  checksum_line <- checksum_lines[grepl(paste0("[[:space:]]", filename, "$"), checksum_lines)]
  if(length(checksum_line) != 1) {
    print(checksum_lines)
    stop(sprintf("Could not find exactly one SHA-256 checksum for %s", filename))
  }

  prev_sha256 <- stringr::str_split_fixed(checksum_line, "[[:space:]]+", 2)[1]
  print(paste("Got previous SHA-256 value to compare against:", prev_sha256))

  f <- RCurl::CFILE(file, mode="wb")
  a <- RCurl::curlPerform(url = url, writedata = f@ref, noprogress=FALSE)
  RCurl::close(f)

  print("Computing SHA-256 for downloaded file")
  new_sha256 <- unname(tools::sha256sum(file))
  print(paste("SHA-256 is:", new_sha256))

  if(new_sha256==prev_sha256) {
    print("SHA-256 matches")
  } else {
    if(file.exists(file)){
      file.remove(file)
    }
    stop("SHA-256 does not match the downloaded file")
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
  tor
}

###############################################
#' Build executable and argv for invoking Bascet directly
#'
#' @param bascetInstance Bascet instance
#' @param params Character vector of Bascet arguments
#'
#' @return List with `command` and `args` elements
#' @noRd
buildBascetInvocation <- function(bascetInstance, params = character()) {
  stopifnot(is.bascet.instance(bascetInstance))

  log_file <- Sys.getenv("BASCET_LOGFILE", unset = NA_character_)
  args <- c(
    if(bascetInstance@logToFile && !is.na(log_file) && nzchar(log_file)) paste0("--log-mode=", log_file),
    paste0("--log-level=", bascetInstance@logLevel),
    params
  )

  prepend_parts <- splitPrependLocal(bascetInstance@prependCmd)
  if(length(prepend_parts) == 0) {
    list(command = bascetInstance@bin, args = args)
  } else {
    list(command = prepend_parts[[1]], args = c(prepend_parts[-1], bascetInstance@bin, args))
  }
}
