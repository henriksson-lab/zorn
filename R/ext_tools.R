
################################################################################
################ FASTQC ######################################################
################################################################################


###############################################
#' Run FASTQC on contigs of all cells.
#' This is a thin wrapper around BascetMapCell
#' 
#' @inheritParams template_BascetFunction
#' @param db description
#' @return TODO
#' @export
BascetMapCellFASTQC <- function( 
    bascetRoot, 
    inputName="contigs",
    outputName="fastqc", 
    #includeCells=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascet_instance=GetDefaultBascetInstance()
){
  BascetMapCell(
    bascetRoot=bascetRoot, 
    withfunction="_fastqc", 
    inputName=inputName, 
    outputName=outputName,
    #includeCells=includeCells
    overwrite=overwrite,
    runner=runner,
    bascet_instance=bascet_instance
  )
}



################################################################################
################ Abricate ######################################################
################################################################################



###############################################
#' Callback function for aggregating ABRicate data for each cell.
#' To be called from BascetAggregateMap
#' 
#' @return TODO
#' @export
aggr.abricate <- function(bascetFile, cellID, bascet_instance){
  tmp <- BascetReadFile(bascetFile, cellID, "abricate.tsv", as="text", bascet_instance=bascet_instance)
  
  #tmp <- readLines("/home/mahogny/github/zorn/test_aggr/abricate/stupid.tsv")
  #tmp <- readLines("/home/mahogny/github/zorn/test_aggr/abricate/salmonella_SRR33219394_ncbi.tsv")
  zz <- textConnection(tmp)
  dat <- read.delim(zz)
  close(zz)
  
  if(nrow(dat)>0){
    dat$cellID <- cellID  #needed? could make this general
  }
  
  #dat <- stringr::str_split_fixed(dat,"\t",15)
  #colnames(dat) <- c("FILE","SEQUENCE","START","END","STRAND","GENE","COVERAGE","COVERAGE_MAP","GAPS","PERC_COVERAGE","PERC_IDENTITY","DATABASE","ACCESSION","PRODUCT","RESISTANCE")
  #dat <- as.data.frame(dat)
  #dat$cellID <- c("A","A","A","B","B","B","B")
  dat
}


###############################################
#' Run Abricate on contigs of all cells.
#' This is a thin wrapper around BascetMapCell
#' 
#' @inheritParams template_BascetFunction
#' @param db description
#' @return TODO
#' @export
BascetMapCellAbricate <- function( 
    bascetRoot, 
    inputName="contigs",
    outputName="abricate", 
    db="ncbi",
    #includeCells=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascet_instance=GetDefaultBascetInstance()
){
  BascetMapCell(
    bascetRoot=bascetRoot, 
    withfunction="_abricate", 
    inputName=inputName, 
    outputName=outputName,
    #includeCells=includeCells
    overwrite=overwrite,
    args = list(DATABASE_DIR=db),
    runner=runner,
    bascet_instance=bascet_instance
  )
}



###############################################
#' List installed databases available for Abricate
#' 
#' @return List of database names
#' @export
ListDatabaseAbricate <- function(
    dbdir,
    bascet_instance=GetDefaultBascetInstance()
) {
  ret <- system(
    paste(
      bascet_instance@prepend_cmd,
      "abricate --list"
    ),
    intern = TRUE
  )
  ret
}




###############################################
#' Aggregate data from Abricate
#' This is a thin wrapper around BascetAggregateMap
#' 
#' @export
BascetAggregateAbricate <- function( 
    bascetRoot, 
    inputName="abricate",
    #cacheFile=NULL, #option
    include_cells=NULL,
    runner=GetDefaultBascetRunner(),
    bascet_instance=GetDefaultBascetInstance()
){
  CountDataFrameToSparseMatrix(MapCellMultiListAsDataFrame(BascetAggregateMap(
    bascetRoot,
    inputName,
    aggr.abricate,
    include_cells=include_cells
  )), "cellID","GENE")
}


################################################################################
################ Bakta #########################################################
################################################################################


###############################################
#' Run Bakta on contigs of all cells.
#' This is a thin wrapper around BascetMapCell
#' 
#' @inheritParams template_BascetFunction
#' @param db description
#' @return TODO
#' @export
BascetMapCellBakta <- function( 
    bascetRoot, 
    inputName="contigs",
    outputName="bakta", 
    db,
    #includeCells=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascet_instance=GetDefaultBascetInstance()
){
  BascetMapCell(
    bascetRoot=bascetRoot, 
    withfunction="_bakta", 
    inputName=inputName, 
    outputName=outputName,
    #includeCells=includeCells
    overwrite=overwrite,
    args = list(DATABASE_DIR=db),
    runner=runner,
    bascet_instance=bascet_instance
  )
}




###############################################
#' Download a database for Bakta
#' @export
DownloadDatabaseBakta <- function(
    dbdir,
    dbtype=c("light","full"),  #todo look up how to handle documentation for this
    bascet_instance=GetDefaultBascetInstance()
) {
  if(file.exists(dbdir)) {
    print("Database already exists; skipping")
  } else {
    system(
      paste(
        bascet_instance@prepend_cmd,
        "bakta_db download --output",
        dbdir,
        "--type",dbtype
      )
    )
  }
}



################################################################################
################ Ariba #########################################################
################################################################################


###############################################
#' Run Ariba on reads of all cells.
#' This is a thin wrapper around BascetMapCell
#' 
#' @inheritParams template_BascetFunction
#' @param db description
#' @return TODO
#' @export
BascetMapCellAriba <- function( 
    bascetRoot, 
    inputName="filtered",
    outputName="ariba", 
    db,
    #includeCells=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascet_instance=GetDefaultBascetInstance()
){
  BascetMapCell(
    bascetRoot=bascetRoot, 
    withfunction="_ariba", 
    inputName=inputName, 
    outputName=outputName,
    #includeCells=includeCells
    overwrite=overwrite,
    args = list(DATABASE_DIR=db),
    runner=runner,
    bascet_instance=bascet_instance
  )
}





###############################################
#' Download database for Ariba
#' 
#' @export
DownloadDatabaseAriba <- function(
    dbdir,
    ref=c("argannot", "card", "ncbi", "megares", "plasmidfinder", "resfinder", "srst2_argannot", "vfdb_core", "vfdb_full", "virulencefinder"),
    bascet_instance=GetDefaultBascetInstance()
) {
  
  tmp <- tempfile()
  
  if(file.exists(dbdir)) {
    print("Database already exists; skipping")
  } else {
    
    ### Get ref
    system(
      paste(
        bascet_instance@prepend_cmd,
        "ariba getref ",ref,tmp,
        "-m", paste0(tmp,".tsv"),
        dbdir
      )
    )
    
    ### prepare the database
    
  }
  
  #todo if above is heavy, make it a slurm job
}

# ariba getref ncbi out.ncbi
# ariba prepareref -f out.ncbi.fa -m out.ncbi.tsv out.ncbi.prepareref

# TODO aggregate scripts for outputs


################################################################################
################ AMRfinder #####################################################
################################################################################



###############################################
#' Run AMRfinder on contigs of all cells.
#' This is a thin wrapper around BascetMapCell
#' 
#' @export
#' @inheritParams template_BascetFunction
#' @param db description
#' @return TODO
#' @export
BascetMapCellAMRfinder <- function( 
    bascetRoot, 
    inputName="contigs",
    outputName="AMRfinder", 
    db,
    #includeCells=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascet_instance=GetDefaultBascetInstance()
){
  BascetMapCell(
    bascetRoot=bascetRoot, 
    withfunction="_amrfinder", 
    inputName=inputName, 
    outputName=outputName,
    #includeCells=includeCells
    overwrite=overwrite,
    args = list(DATABASE_DIR=db),
    runner=runner,
    bascet_instance=bascet_instance
  )
}



###############################################
#' Download a database for AMRfinder
#' 
#' @export
DownloadDatabaseAMRfinder <- function(
    dbdir,
    bascet_instance=GetDefaultBascetInstance()
) {
  if(file.exists(dbdir)) {
    print("Database already exists; skipping")
  } else {
    system(
      paste(
        bascet_instance@prepend_cmd,
        "amrfinder_update -d",
        dbdir
      )
    )
  }
}




