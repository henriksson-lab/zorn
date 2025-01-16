


################################################################################
################ Reading of Bascet files #######################################
################################################################################


#' A bascet, along with all the shards
#' @export
setClass("Bascet", slots=list(
  num_shards="numeric",
  files="character",
  cellmeta="ANY"
)
) 


###############################################
#' Get list of cells in a bascet
#' @return TODO
#' @export
BascetCellNames <- function(bascetRoot, bascetName){
  
  #Can support BAM later as well
  shards <- detect_shards_for_file(bascetRoot,bascetName)
  if(length(shards)==0){
    stop(paste("Could not find file ",bascetName))
  }
  
  cellnames <- list()
  for(i in seq_along(shards)){
    cur_file <- file.path(bascetRoot,shards[i])
    if(stringr::str_ends(cur_file, fixed(".zip"))) {
      allfiles <- unzip(cur_file,list = TRUE)$Name
    } else if(stringr::str_ends(cur_file, fixed(".tirp.gz"))) {
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



system(paste("tabix","--list-chroms", "/home/mahogny/github/bascet/testdata/filtered.0.tirp.gz"),intern=TRUE)


if(FALSE){
  BascetCellNames("/home/mahogny/jupyter/zorn/test/data", "shard")
}



#' Open a Bascet, prepare it for reading.
#' 
#' The current code is based on pure R, but more efficient calls can be made
#' in the future. We thus advise against direct zip-file manipulation and
#' do not guarantee future support for this
#' @return TODO
#' @export
OpenBascet <- function(bascetRoot, bascetName){
  
  shards <- detect_shards_for_file(bascetRoot,bascetName)
  num_shards <- length(shards)
  cellname_coord <- BascetCellNames(bascetRoot, bascetName)
  
  new("Bascet", 
      num_shards=num_shards, 
      files=file.path(bascetRoot, shards), 
      cellmeta=BascetCellNames(bascetRoot, bascetName))
}




#' Read one file from a Bascet
#' 
#' This can be made faster by, e.g., once and for all reading the location of
#' all objects in the file
#' 
#' @return TODO
#' @export
BascetReadFile <- function(bascetFile, cellID, filename, as=c("tempfile"), bascet_instance=bascet_instance.default, verbose=FALSE){
  
  ## Check if the cell is present at all
  cellmeta <- bascetFile@cellmeta[bascetFile@cellmeta$cell==cellID,,drop=FALSE]
  if(nrow(cellmeta)==0){
    stop("Cell not present in file")
  }
  
  if(as=="pipe"){
    #bascet pipe <<<file  <<<cellname
    stop("not implemented")    
  } else if(as=="tempfile"){
    
    #Need a directory to unzip to
    tname.dir <- tempfile()
    dir.create(tname.dir)
    
    #Extract this zip file and then check that it worked
    name_of_zip <- bascetFile@files[cellmeta$shard+1]
    extract_files <- file.path(cellID,filename)
    
    if(verbose){
      print(name_of_zip)
      print(extract_files)
      print(tname.dir)
    }
    
    
    if(FALSE){
      ########### Pure R version. does not support our zip format
      unzip(name_of_zip, files=extract_files, exdir=tname.dir) ### 666 cannot operate on our rust files. 
      
      name_of_outfile <- file.path(tname.dir, cellID, filename) ## should be checked
      if(!file.exists(name_of_outfile)){
        stop(paste("unzip failed to produce expected ",name_of_outfile))
      }
      
      #Move the output file to a new temp file, such that the callback function need not delete whole directory
      tname.out <- tempfile()
      file.rename(name_of_outfile, tname.out)
      
      #Recursive delete. to be safe, remove expected files
      file.remove(file.path(tname.dir, cellID))
      file.remove(tname.dir)
      
    } else {
      ########### Use Bascet to unzip
      #cargo +nightly run extract -i /Users/mahogny/Desktop/rust/hack_robert/testdata/quast.zip  -o /Users/mahogny/Desktop/rust/hack_robert/testdata/out.temp -b a  -f report.txt
      tname.out <- tempfile()
      
      cmd = paste(
        bascet_instance@bin, 
        "extract -i",name_of_zip, 
        #"-t", bascet_instance@tempdir,
        "-o",tname.out,
        "-b",cellID,
        "-f",filename
      )
      system(cmd)#, show.output.on.console = TRUE)
      #      files_in[$TASK_ID] --o files_out[$TASK_ID] -f ", withfunction),
    }
    
    
    #Return temp file location
    return(tname.out)
  } else {
    stop("Unsupported output format")
  }
}










#' Read one file from a Bascet
#' 
#' This can be made faster by, e.g., once and for all reading the location of
#' all objects in the file
#' 
#' @return TODO
#' @export
BascetListFilesForCell <- function(bascetFile, cellID, bascet_instance=bascet_instance.default){
  
  ## Check if the cell is present at all
  cellmeta <- bascetFile@cellmeta[bascetFile@cellmeta$cell==cellID,,drop=FALSE]
  if(nrow(cellmeta)==0){
    stop("Cell not present in file")
  }

  #Extract this zip file and then check that it worked
  name_of_zip <- bascetFile@files[cellmeta$shard+1]

  res <- unzip(name_of_zip, list=TRUE)
  div <- stringr::str_split_fixed(res$Name,stringr::fixed("/"),2)
  data.frame(
    cell=div[,1],
    file=div[,2],
    size=res$Length
  )
}






################################################################################
################ Histogram functions ###########################################
################################################################################


#' @return TODO
#' @export
ReadHistogram <- function(bascetRoot, inputName, bascet_instance=bascet_instance.default){
  
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


#' @return TODO
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
################ Atrandi-specific #######################################
################################################################################


#' @return TODO
#' @export
AtrandiBarcodeStats <- function(
    bascetRoot, 
    inputName="debarcoded", 
    bascet_instance=bascet_instance.default){

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


#' Helper function: Figure out which shards belong together given root input name and extension
#' i.e. root/name.##.ext
detect_shards_for_file <- function(bascetRoot, inputName){
  allfiles <- list.files(bascetRoot)
  allzip <- na.omit(stringr::str_match(allfiles, paste0(inputName,"\\.[0123456789]+\\.","zip")))
  alltirp <- na.omit(stringr::str_match(allfiles, paste0(inputName,"\\.[0123456789]+\\.","tirp.gz")))
  unique(c(allzip,alltirp)) #tirp.gz can sneak in
}

#detect_shards_for_file("~/jupyter/zorn/test","fakein","zip")

#' Helper function: Generate suitable output filenames according to shard system
#' i.e. root/name.##.ext
make_output_shard_names <- function(bascetRoot, outputName, ext, num_shards){
  file.path(bascetRoot, paste0(outputName,".",seq_len(num_shards),".",ext))
}

#make_output_shard_names("/","bar","zip", 5)

