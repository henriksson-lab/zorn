
################################################################################
################ FASTQC ########################################################
################################################################################


###############################################
#' Run FASTQC on reads of all cells.
#' This is a thin wrapper around BascetMapCell
#' 
#' @inheritParams template_BascetFunction
#' @param db description
#' @return TODO
#' @export
BascetMapCellFASTQC <- function( 
    bascetRoot, 
    inputName="filtered",
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




###############################################
#' Aggregate data from FASTQC
#' This is a thin wrapper around BascetAggregateMap
#' 
#' @export
BascetAggregateFASTQC <- function( 
    bascetRoot, 
    inputName="fastqc",
    #cacheFile=NULL, #option
    include_cells=NULL,
    runner=GetDefaultBascetRunner(),
    bascet_instance=GetDefaultBascetInstance()
){
  CountDataFrameToSparseMatrix(BascetAggregateMap(
    bascetRoot,
    inputName,
    aggr.fastqc,
    include_cells=include_cells
  ))
}



###############################################
#' Callback function for aggregating FASTQC data for each cell.
#' To be called from BascetAggregateMap
#' 
#' @return TODO
#' @export
aggr.fastqc <- function(bascetFile, cellID, bascet_instance){
  
  #tmp <- readLines("/home/mahogny/github/zorn/test_aggr/abricate/stupid.tsv")
  #tmp <- readLines("/home/mahogny/github/zorn/test_aggr/abricate/salmonella_SRR33219394_ncbi.tsv")
  #lines <- readLines("/home/mahogny/github/zorn/test_aggr/fastqc/new/r1_fastqc_data.txt")
  #internal_parse_fastqc_data(lines)
  
  #Handle both R1 and R2
  lines_r1 <- BascetReadFile(bascetFile, cellID, "r1_fastqc_data.txt", as="text", bascet_instance=bascet_instance)
  lines_r2 <- BascetReadFile(bascetFile, cellID, "r2_fastqc_data.txt", as="text", bascet_instance=bascet_instance)

  
  #lines_r1 <- readLines("/home/mahogny/github/zorn/test_aggr/fastqc/new/r1_fastqc_data.txt")
  
  dat_r1 <- internal_parse_fastqc_data(lines_r1)
  dat_r2 <- internal_parse_fastqc_data(lines_r2)
  
  #Return both R1 and R2
  list(
    r1=dat_r1, 
    r2=dat_r2
  )
}

#Parse one FASTQC report file, divide into tables for each section
internal_parse_fastqc_data <- function(lines){
  
  #Separate state of section
  all_section_states <- list()
  
  #Put each separate section in a list
  list_sections <- list()
  lines <- lines[lines!=">>END_MODULE" & lines!=""]
  section_start <- which(stringr::str_detect(lines,stringr::fixed(">>")))
  section_end <- c(section_start[-1], length(lines))
  for(i in 1:length(section_start)) {
    subsection <- lines[section_start[i]:(section_end[i]-1)]
    
    zz <- textConnection(subsection[-1])
    dat <- read.delim(zz)
    close(zz)
    
    
    #Divide section and state
    section_name_state <- stringr::str_sub(subsection[1],3)
    section_name_state <- stringr::str_split_fixed(section_name_state, "\t",2)
    
    #to_store <- section_name_state[,2]
    #names(to_store) <- section_name_state[,1]
    #print(777)
    #print(t(to_store))
    all_section_states[[i]] <- data.frame(state=section_name_state[,1], val=section_name_state[,2])
    #rbind(
    #  all_section_states, 
    #  to_store # section_name_state
    #)
    
    #Store this data frame in current section
    section_name <- section_name_state[1]
    list_sections[[section_name]] <- dat
  }

  all_section_states <- do.call(rbind,all_section_states)
  all_section_states <- data.frame(
    row.names=all_section_states$state,
    val=all_section_states$val
    )
  
  
  all_section_states <- t(all_section_states)
  print(all_section_states)
  
  #Name section-state table and include it as well  
  #colnames(all_section_states) <- c("section","state")
  #all_section_states <- as.data.frame(all_section_states)
  list_sections[["section_state"]] <- all_section_states
  
  # list_sections <- internal_parse_fastqc_data(readLines("/home/mahogny/github/zorn/test_aggr/fastqc/new/r1_fastqc_data.txt"))
  
  #print(666)
  #print(all_section_states)
  
  list_sections
}




###############################################
#' Callback function for aggregating FASTQC data for each cell.
#' To be called from BascetAggregateMap
#' 
#' @param readnum 1 or 2, for R1 and R2
#' @param use_browser Use operating system browser to open file
#' @return TODO
#' @export
ShowFASTQCforCell <- function(
    bascetFile, 
    cellID, 
    readnum="1",
    use_browser=FALSE,
    verbose=FALSE
){
  #Open the bascet file, get the HTML report
  if(verbose){
    print("Creating extract streamer session")
  }
  bascetFile <- OpenBascet(bascetRoot, bascetName)
  if(verbose){
    print("Extract streamer session ok")
  }
  lines <- BascetReadFile(bascetFile, cellID, paste0("r",readnum,"_fastqc.html"), as="text", bascet_instance=bascet_instance)
  #lines <- readLines("/home/mahogny/github/zorn/test_aggr/fastqc/new/r1_fastqc.html")
  CloseBascet(bascetFile)
  
  #RStudio is picky. File to show really has to be created this way
  dir <- tempfile()
  dir.create(dir)
  htmlFile <- file.path(dir, "index.html")
  writeLines(lines, htmlFile)

  #Show the file
  viewer <- getOption("viewer")
  if (!is.null(viewer) & !use_browser) {
    viewer(htmlFile)
  } else {
    utils::browseURL(htmlFile)
  }
  Invisible()
}


###############################################
# Helper for FASTQC: get a data.frame across all cells for a section and read
###############################################
#' Get a data frame for one type of FASTQ statistics across across all cells
#' 
#' @param readnum 1 or 2, for R1 and R2
#' @return TODO
#' @export

GetFASTQCassembledDF <- function(
    aggr_fastq_data, 
    section, 
    readnum
) {
  
  internal_fastqc_getread_in_list <- function(lst,readnum){
    lapply(lst, function(s) s[[paste0("r",readnum)]])
  }

  internal_fastqc_add_cellid_to_list <- function(lst){
    lapply(names(lst), function(x) {
      temp <- as.data.frame(lst[[x]]) #to be on the safe side; fixed one bug
      temp$cellID <- x
      temp
    })
  }
  
  list_oneread <- internal_fastqc_getread_in_list(aggr_fastq_data,readnum) #Get data for given read
  #print(list_oneread)
  
  list_oneread_section <- lapply(list_oneread, function(s) s[[section]]) #Get sections for each cell
  
  df_section <- do.call(rbind, internal_fastqc_add_cellid_to_list(list_oneread_section))
  df_section
}


###############################################
#' From aggregated FASTQC data, plot adapter content
#' 
#' TODO if giving many cells
#' 
#' @param readnum 1 or 2, for R1 or R2
#' @export
PlotFASTQCadapterContent <- function(
    aggr_fastq_data,
    readnum
) {
  df_section <- GetFASTQCassembledDF(aggr_fastq_data,"Adapter Content",readnum)
  xlab_unit <- unique(df_section$X.Position)
  
  df_section$any_adapter <- rowSums(df_section[,!(colnames(df_section) %in% c("X.Position","cellID"))])
  df_section$X.Position <- factor(df_section$X.Position, levels=xlab_unit) 
  
  #Possible subsets; but would need one plot window for each
  #Illumina.Universal.Adapter Illumina.Small.RNA.3..Adapter Illumina.Small.RNA.5..Adapter Nextera.Transposase.Sequence SOLID.Small.RNA.Adapter cellID
  
  ggplot(df_section, aes(X.Position, any_adapter, group=cellID, color=cellID)) + 
    geom_line() + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylab("% Adapters")  
}




###############################################
#' From aggregated FASTQC data, get overall statistics for overlay on UMAP etc
#' 
#' TODO argument if giving many cells
#' 
#' @param readnum 1 or 2, for R1 or R2
#' @export
GetFASTQCpassfailStats <- function(
    aggr_fastq_data,
    readnum
) {
  df <- GetFASTQCassembledDF(aggr_fastq_data,"section_state",readnum)
  rownames(df) <- df$cellID
  df <- df[colnames(df)!="cellID",drop=FALSE]  
  #Could do both r1 and r2
  df  
}





if(FALSE){
  list_sections <- internal_parse_fastqc_data(readLines("/home/mahogny/github/zorn/test_aggr/fastqc/new/r1_fastqc_data.txt"))
  
  testpair <- list(
    r1=list_sections,
    r2=list_sections
  )
  
  testlist <- list(
    "A"=testpair,
    "B"=testpair
  )
  
  # [1] "Basic Statistics"
  # [1] "Per base sequence quality"
  # [1] "Per sequence quality scores"
  # [1] "Per base sequence content"
  # [1] "Per sequence GC content"
  # [1] "Per base N content"
  # [1] "Sequence Length Distribution"
  # [1] "Sequence Duplication Levels"
  # [1] "Overrepresented sequences"
  # [1] "Adapter Content"
  
  # PlotFASTQCadapterContent(testlist,1)
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
    dat$cellID <- cellID  #needed? could make this general  ---- think no longer used!
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
    ref="ncbi", #c("argannot", "card", "ncbi", "megares", "plasmidfinder", "resfinder", "srst2_argannot", "vfdb_core", "vfdb_full", "virulencefinder"),
    bascet_instance=GetDefaultBascetInstance()
) {
  
  tmp <- tempfile()
  tmp.fa <- paste0(tmp,".fa")
  tmp.tsv <- paste0(tmp,".tsv")
  out.prepareref <- file.path(dbdir,"out.prepareref")  #out.ncbi.prepareref
  
  #dir.exists()  better?
  
  if(file.exists(dbdir)) {
    print("Database already exists; skipping")
  } else {
    dir.create(dbdir, recursive=TRUE)
    
    # ariba getref ncbi out.ncbi
    # ariba prepareref -f out.ncbi.fa -m out.ncbi.tsv out.ncbi.prepareref
    
    ### Get ref
    system(
      paste(
        bascet_instance@prepend_cmd,
        "ariba getref ",ref, tmp, 
        dbdir
      )
    )
    
    ### prepare the database
    system(
      paste(
        bascet_instance@prepend_cmd,
        "ariba prepareref -f ",tmp.fa,
        "-m", tmp.tsv, out.prepareref
        #        dbdir
      )
    )
    
    stop("todo implement")
  }
  
  #todo if above is heavy, make it a slurm job
}


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


###




###############################################
#' Callback function for aggregating ABRicate data for each cell.
#' To be called from BascetAggregateMap
#' 
#' @return TODO
#' @export
aggr.amrfinder <- function(bascetFile, cellID, bascet_instance){
  tmp <- BascetReadFile(bascetFile, cellID, "amrfinder.tsv", as="text", bascet_instance=bascet_instance)
  
  #tmp <- readLines("/home/mahogny/github/zorn/test_aggr/amrfinder/salmonella_SRR33219394_amrfinder.tsv")
  zz <- textConnection(tmp)
  dat <- read.delim(zz)
  close(zz)
  
  #if(nrow(dat)>0){
  #  dat$cellID <- cellID  #needed? could make this general  -- think no longer used
  #}
  
  dat
}


"Element.symbol"


###############################################
#' Aggregate data from AMRfinder
#' This is a thin wrapper around BascetAggregateMap
#' 
#' @export
BascetAggregateAMRfinder <- function( 
    bascetRoot, 
    inputName="AMRfinder",
    #cacheFile=NULL, #option
    include_cells=NULL,
    get_column="Element.symbol",  #Element.name is an option, if one also want full name
    runner=GetDefaultBascetRunner(),
    bascet_instance=GetDefaultBascetInstance()
){
  CountDataFrameToSparseMatrix(MapCellMultiListAsDataFrame(BascetAggregateMap(
    bascetRoot,
    inputName,
    aggr.abricate,
    include_cells=include_cells
  )), "cellID",get_column)
}

