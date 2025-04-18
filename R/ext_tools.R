


###############################################
#' ....
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
#' ....
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
#' ....
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



###############################################
#' ....
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






###############################################
#' ....
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



