


################################################################################
################ Reading of Bascet files #######################################
################################################################################


###############################################
#' A bascet, along with all the shards
#' @export
setClass("Bascet", slots=list(
  num_shards="numeric",
  files="character",
  cellmeta="ANY",
  streamer="ANY"
)
) 


###############################################
#' Get list of cells in a Bascet 
#' 
#' @inheritParams template_BascetFunction
#' @param bascetName Name of the bascet
#' @return Vector of cell names as strings
#' @export
BascetCellNames <- function(
    bascetRoot, 
    bascetName
){
  
  #Can support BAM later as well
  shards <- detect_shards_for_file(bascetRoot,bascetName)
  if(length(shards)==0){
    stop(paste("Could not find file ",bascetName))
  }
  
  cellnames <- list()
  for(i in seq_along(shards)){
    cur_file <- file.path(bascetRoot,shards[i])
    if(stringr::str_ends(cur_file, stringr::fixed(".zip"))) {
      allfiles <- unzip(cur_file,list = TRUE)$Name
    } else if(stringr::str_ends(cur_file, stringr::fixed(".tirp.gz"))) {
      allfiles <- system(paste("tabix","--list-chroms", cur_file),intern=TRUE)
    } else {
      stop(paste("Cannot list cells for"),cur_file)
    }
    df <- data.frame(cell=unique(stringr::str_split_i(allfiles,"/",1)))   
    df$shard <- i - 1
    cellnames[[i]] <- df
  }
  do.call(rbind, cellnames)
}





if(FALSE){
  system(paste("tabix","--list-chroms", "/home/mahogny/github/bascet/testdata/filtered.0.tirp.gz"),intern=TRUE)
  BascetCellNames("/home/mahogny/jupyter/zorn/test/data", "shard")
}



###############################################
#' Open a Bascet, prepare it for reading individual files
#' 
#' 
#' TODO The current code is based on pure R, but more efficient calls can be made
#' in the future. We thus advise against direct zip-file manipulation and
#' do not guarantee future support for this
#' 
#' @inheritParams template_BascetFunction
#' @return A handle to a Bascet
#' @export
OpenBascet <- function(
    bascetRoot, 
    bascetName,
    bascet_instance=GetDefaultBascetInstance()
){
  
  shards <- detect_shards_for_file(bascetRoot,bascetName)
  num_shards <- length(shards)
  cellname_coord <- BascetCellNames(bascetRoot, bascetName)
  
  streamer <- extractstreamer_start(bascet_instance = bascet_instance)
  
  new("Bascet", 
      num_shards=num_shards, 
      files=file.path(bascetRoot, shards), 
      cellmeta=BascetCellNames(bascetRoot, bascetName),
      streamer=streamer
  )
}





###############################################
#' Close a Bascet file.
#' 
#' This should always be performed to avoid memory leaks
#' 
#' @inheritParams template_BascetFunction
#' @return A handle to a Bascet
#' @export
CloseBascet <- function(
    bascetFile
){
  extractstreamer_exit(bascetFile@streamer)
  invisible()
}





###############################################
#' Read one file from a Bascet
#' 
#' TODO This can be made faster by, e.g., once and for all reading the location of
#' all objects in the file
#' 
#' @inheritParams template_BascetFunction
#' @param filename description
#' @param as description
#' @return Name of a temporary file where the read content is stored
#' @export
BascetReadFile <- function(
    bascetFile, 
    cellID, 
    filename, 
    as=c("tempfile"), 
    bascet_instance=GetDefaultBascetInstance(), 
    verbose=FALSE
){
  
  ## Check if the cell is present at all
  cellmeta <- bascetFile@cellmeta[bascetFile@cellmeta$cell==cellID,,drop=FALSE]
  if(nrow(cellmeta)==0){
    stop("Cell not present in file")
  }
  
  
  name_of_zip <- bascetFile@files[cellmeta$shard+1]
  extract_files <- file.path(cellID,filename)
  
  
  if(as=="text"){
    
    if(extractstreamer_open(bascetFile@streamer, name_of_zip)==0) {
      ret <- extractstreamer_showtext(bascetFile@streamer, extract_files)
      return(ret)
    } else {
      print(paste("Could not open file",name_of_zip))
      return(NULL)
    }  #is this ok? check return value TODO
    
  } else if(as=="tempfile"){
    
    #Need a directory to unzip to
    tname.dir <- tempfile()
    dir.create(tname.dir)
    
    #Extract this zip file and then check that it worked
    if(verbose){
      print(name_of_zip)
      print(extract_files)
      print(tname.dir)
    }
    
    
    ### TODO use new API
    
    
    ########### Use Bascet to unzip
    #cargo +nightly run extract -i /Users/mahogny/Desktop/rust/hack_robert/testdata/quast.zip  -o /Users/mahogny/Desktop/rust/hack_robert/testdata/out.temp -b a  -f report.txt
    
    ret <- extractstreamer_open(bascetFile@streamer, name_of_zip)
    if(ret==0) {
      tname.out <- tempfile()
      ret <- extractstreamer_extract_to(bascetFile@streamer, extract_files, tname.out)
      if(ret==0) {
        #Return temp file location
        return(tname.out)
      } else {
        print(paste("cell", cellID,": Failed to get file", filename, " from shard", name_of_zip, ":", ret))
        return(NULL)      
      }
    } else {
      print(paste("cell", cellID,": Failed to open shard", name_of_zip, ":", ret))
      return(NULL)      
    }

  } else {
    stop("Unsupported output format")
  }
}










###############################################
#' List files for a cell in a Bascet
#' 
#' This can be made faster by, e.g., once and for all reading the location of
#' all objects in the file
#' 
#' @param bascetFile Bascet file object
#' @param cellID Name of the cell
#' @param bascet_instance A Bascet instance
#' @return A data.frame with list of all the files
#' @export
BascetListFilesForCell <- function(
    bascetFile, 
    cellID, 
    bascet_instance=GetDefaultBascetInstance(),
    super_verbose=FALSE
){
  
  ## Check if the cell is present at all
  cellmeta <- bascetFile@cellmeta[bascetFile@cellmeta$cell==cellID,,drop=FALSE]
  if(nrow(cellmeta)==0){
    stop("Cell not present in file")
  }

  #Extract this zip file and then check that it worked
  name_of_zip <- bascetFile@files[cellmeta$shard+1]

  ret <- extractstreamer_open(bascetFile@streamer, name_of_zip, super_verbose)  #is this ok? check return value TODO
  if(ret==0){
    res <- extractstreamer_ls(bascetFile@streamer)
    
    div <- stringr::str_split_fixed(res$Name,stringr::fixed("/"),2)
    data.frame(
      cell=div[,1],
      file=div[,2],
      size=res$Length
    )    
  } else {
    print("Could not open shard",name_of_zip)
    data.frame()
  }
  

}






################################################################################
################ Histogram functions ###########################################
################################################################################


###############################################
#' 
#' Read the count histogram associated with a Bascet.
#' Not all Bascets have one, but it is typically produced after debarcoding
#' 
#' @param bascet_instance A Bascet instance
#' @return Histogram as a data.frame
#' @export
ReadHistogram <- function(
    bascetRoot, 
    inputName, 
    bascet_instance=GetDefaultBascetInstance()
){
  
  #Get all the TIRPs, sum up the reads  
  inputFiles <- detect_shards_for_file(bascetRoot, inputName)
  print(inputFiles)
  
  list_hist <- list()
  for(f in inputFiles) {
    
    hist_f <- paste0(file.path(bascetRoot, f),".hist") ## only support TIRP for now
    dat <- read.csv(hist_f, sep="\t")
    list_hist[[f]] <- dat
  }
  dat <- do.call(rbind, list_hist)
  colnames(dat) <- c("cellid","count")
  
  dat <- sqldf::sqldf("select cellid, sum(count) as count from dat group by cellid")
  dat <- dat[order(dat$count, decreasing=TRUE),]
  dat 
}


###############################################
#' Plot a histogram, loaded by ReadHistogram
#' 
#' @return A ggplot object
#' @export
PlotHistogram <- function(dat){
  dat <- dat[order(dat$count, decreasing=TRUE),]
  dat$index <- 1:nrow(dat)
  
  ggplot2::ggplot(dat, ggplot2::aes(index,count)) + 
    ggplot2::geom_line() +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::theme_bw()+
    ggplot2::xlab("Cell index") +
    ggplot2::ylab("Number of reads/cell")
  
}

if(FALSE){
  
  head(ReadHistogram("/home/mahogny/github/bascet/testdata","out_complete"))
  PlotHistogram(ReadHistogram("/home/mahogny/github/bascet/testdata","out_complete"))
  
}






################################################################################
################ Atrandi-specific ##############################################
################################################################################


###############################################
#' Given a Bascet, produce a matrix showing for each combinatorial barcode,
#' how many times it occurs across the cells. Presented as a 96-well plate matrix
#' 
#' This command assumes that cells are named as follows: well1_well2_well3_well4,
#' where e.g. well1 is in format G12
#' 
#' @return TODO
#' @export
AtrandiBarcodeStats <- function(
    bascetRoot, 
    inputName="debarcoded", 
    bascet_instance=GetDefaultBascetInstance()
){

  #Get frequency of each barcode  
  cb <- as.data.frame(stringr::str_split_fixed(BascetCellNames(bascetRoot, inputName)$cell,"_",4))
  df <- as.data.frame(table(unlist(cb)))
  colnames(df) <- c("well","cnt")
  
  #Split up into row and columns. Present as a matrix
  df$col <- stringr::str_sub(df$well,1,1)
  df$row <- stringr::str_sub(df$well,2)
  mat <- reshape2::acast(df, col ~ row, value.var = "cnt")
  mat[order(rownames(mat)), order(as.integer(colnames(mat)))]
}



################################################################################
################ helper functions ##############################################
################################################################################


###############################################
#' Helper function: Figure out which shards belong together given root input name and extension
#' i.e. root/name.##.ext
detect_shards_for_file <- function(
    bascetRoot, 
    inputName
){
  allfiles <- list.files(bascetRoot)

  allzip <- stringr::str_detect(allfiles, paste0("^", inputName,"\\.[0123456789]+\\.","zip", "$"))
  alltirp <- stringr::str_detect(allfiles, paste0("^", inputName,"\\.[0123456789]+\\.","tirp\\.gz", "$"))
  allbam <- stringr::str_detect(allfiles, paste0("^", inputName,"\\.[0123456789]+\\.","bam", "$"))
  allcram <- stringr::str_detect(allfiles, paste0("^", inputName,"\\.[0123456789]+\\.","cram", "$"))
  all_kraken <- stringr::str_detect(allfiles, paste0("^", inputName,"\\.[0123456789]+\\.","kraken_out", "$"))
  
  #TODO: rename to kraken5?
  all_kraken_counts <- 
    stringr::str_detect(allfiles, paste0("^", inputName,"\\.[0123456789]+\\.","counts\\.hdf5", "$")) | #remove in future
    stringr::str_detect(allfiles, paste0("^", inputName,"\\.[0123456789]+\\.","kraken5", "$"))  

  all_feature_counts <- 
    stringr::str_detect(allfiles, paste0("^", inputName,"\\.[0123456789]+\\.","hd5", "$")) | #todo remove
    stringr::str_detect(allfiles, paste0("^", inputName,"\\.[0123456789]+\\.","h5", "$"))  # better name?

  
  allfq <- stringr::str_detect(allfiles, paste0("^", inputName,"\\.[0123456789]+\\.","fq\\.gz", "$"))
  allfq_r1 <- stringr::str_detect(allfiles, paste0("^", inputName,"\\.[0123456789]+\\.","R1\\.fq.gz", "$"))

  allfiles[
    allzip |
    alltirp |
    allbam |
    allcram |
    allfq |
    allfq_r1 |
    all_kraken_counts |
    all_feature_counts |
    all_kraken
  ]
}


###############################################
#' Helper function: Generate suitable output filenames according to shard system
#' i.e. root/name.##.ext
make_output_shard_names <- function(
    bascetRoot, 
    outputName, 
    ext, 
    num_shards
){
  file.path(bascetRoot, paste0(outputName,".",seq_len(num_shards),".",ext))
}


if(FALSE){
  #detect_shards_for_file("~/jupyter/zorn/test","fakein","zip")
  #make_output_shard_names("/","bar","zip", 5)
}



























################################################################################
################ Extract Streamer API ##########################################
################################################################################



###############################################
#' extract streamer: create an instance
#' 
#' @return The process
extractstreamer_start <- function(
    fname=NULL,
    verbose=FALSE,
    bascet_instance=GetDefaultBascetInstance()
) {
  
  #Assemble the command
  all_cmd <- stringr::str_trim(paste(
    bascet_instance@prepend_cmd,
    bascet_instance@bin,
    "extract-stream",
    if(!is.null(fname)) c("-i",fname)
  ))
  
  #Get home variable. will be "" if empty
  env_home <- Sys.getenv("HOME") 

  #Somehow ~ did not get expanded... so do the work if needed
  all_cmd <- stringr::str_replace_all(all_cmd, stringr::fixed("~"), env_home)
  
  #processx wants the command delivered one argument at a time
  all_cmd_split <- stringr::str_split(all_cmd, " ")[[1]]
  all_cmd_split <- all_cmd_split[all_cmd_split!=""] #not sure if needed; but helped make it run
  if(verbose){
    print(all_cmd_split)
  }
  p <- processx::process$new(
    all_cmd_split[1], 
    all_cmd_split[-1],
    stdin="|", stdout = "|", stderr = "|"
  )
  all_out <- c()
  while(TRUE){
    if(!p$is_alive()){
      print(p$read_all_error())
      stop("Streamer unexpectedly died")
    }
    
    newlines <- p$read_output_lines()
    all_out <- c(all_out, newlines)
    if(verbose){
      print(newlines)
    }
    
    if(length(all_out)>0) {
      last_line <- all_out[length(all_out)] 
      if(last_line=="ready"){
        return(p)
      } else if(stringr::str_starts(last_line,"error")) {
        print(paste("error from extract streamer start:",last_line))
        break
      }
    }
  }
  p
}


###############################################
#' extract streamer: list all files in current zip-file
#' 
#' @return list of files
extractstreamer_ls <- function(p, super_verbose=FALSE){
  p$write_input("ls\n")
  #Figure out how many lines to get
  newlines <- extractstreamer_read_one_line(p)
  n_lines <- as.integer(newlines)
  extractstreamer_read_n_lines(p, n_lines, super_verbose)
}





###############################################
#' extract streamer: end an instance.
#' the object should no longer be used after calling this function
#' 
extractstreamer_exit <- function(p){
  p$write_input("exit\n")
}


###############################################
#' extract streamer: helper function to read one line
#' 
#' @return the line
extractstreamer_read_one_line <- function(p,super_verbose){
  while(TRUE){
    newlines <- p$read_output_lines(n = 1)
    if(length(newlines)>0) {
      print(paste("superverbose, got line: ",newlines))
      return(newlines)
    }
  }
}


###############################################
#' extract streamer: helper function to read N lines
#' 
#' @return all the lines
extractstreamer_read_n_lines <- function(p, n_lines, super_verbose=FALSE){
  all_out <- c()
  while(length(all_out) < n_lines){
    if(!p$is_alive()){
      print(p$read_all_error())
      stop("process unexpectedly died")
    }
    newlines <- p$read_output_lines()
    if(super_verbose){
      print("got more lines:")
      print(newlines)
      print(paste("len all out", length(all_out)))
      print(paste("n_lines", n_lines))
    }
    all_out <- c(all_out, newlines)
  }
  if(super_verbose){
    print("return from extractstreamer_read_n_lines")
  }
  all_out
}

###############################################
#' extract streamer: get content of file, assumed to be text (or this function crashes)
#' 
#' @return The text, divided by line
extractstreamer_showtext <- function(p, fname, super_verbose) {
  p$write_input(paste0("showtext ",fname,"\n"))

  #Figure out how many lines to get, then get them
  newlines <- extractstreamer_read_one_line(p, super_verbose)
  if(stringr::str_starts(newlines, "error")){
    #print(newlines)
    #Return nothing
    NULL
  } else {
    #Get all the lines
    n_lines <- as.integer(newlines)
    extractstreamer_read_n_lines(p, n_lines, super_verbose)
  }
}


###############################################
#' extract streamer: set which file is open
#' 
#' @return 0 if ok
extractstreamer_open <- function(
    p, 
    fname, 
    super_verbose=FALSE
) {
  p$write_input(paste0("open ",fname,"\n"))
  #Figure out how many lines to get, then get them
  newlines <- extractstreamer_read_one_line(p,super_verbose)
  if(newlines=="ok"){
    0
  } else {
    print(paste("from extract streamer:",newlines))
    1
  }
}


###############################################
#' extract streamer: extract to external file
#' @return 0 if ok
extractstreamer_extract_to <- function(
    p, 
    fname, 
    outpath
) {
  p$write_input(paste0("extract_to ",fname," ",outpath,"\n"))
  newlines <- extractstreamer_read_one_line(p) #done or error
  if(newlines=="ok"){
    0
  } else {
    print(newlines)
    1
  }
}




################################################################################
########### Testing ############################################################
################################################################################


if(FALSE){
  system("singularity run /home/mahogny/github/bascet/singularity/bascet.sif ls /home/mahogny/github/bascet/testdata/minhash.0.zip")
  system("singularity run /home/mahogny/github/bascet/singularity/bascet.sif ls /home/mahogny/foo/minhash.0.zip")
  system("ls /home/mahogny/github/bascet/testdata/minhash.0.zip")
  system("singularity run /home/mahogny/github/bascet/singularity/bascet.sif bascet extract-stream -i /home/mahogny/github/bascet/testdata/minhash.0.zip")
}


if(FALSE){
  #p <- extractstreamer_start("/home/mahogny/github/bascet/testdata/minhash.0.zip")
  p <- extractstreamer_start()
  extractstreamer_open(p, "/home/mahogny/foo/minhash.0.zip")
  extractstreamer_ls(p)
  extractstreamer_showtext(p,"C1_B4_B9_A12/cellmap.log")
  for(i in 1:100){
    extractstreamer_showtext(p,"C1_B4_B9_A12/minhash.txt")
  }
  
  extractstreamer_extract_to(p,"C1_B4_B9_A12/minhash.txt","/home/mahogny/temp.bar")
  
}






