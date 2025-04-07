

###############################################
#' Function template, where basic parameter documentation can be obtained from
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascet_instance Configuration for how to run the Bascet Rust API
#' @param num_output_shards Number of output shards, i.e., how many output files to split the data into? (>=1)
#' @param includeCells List of cells to process
#' @param bascetFile Handle for an opened Bascet file
#' @param cellID ID of the cell (string)
#' @param verbose description
#' 
#' @return A job to be executed
template_BascetFunction <- function(
    bascetRoot, 
    runner, 
    bascet_instance=bascet_instance.default){}




###############################################
#' Detect metadata for raw input FASTQ files
#' 
#' _R1 -- common from illumina sequencer
#' SRR****_1.fastq.gz -- typical from SRA
#' 
#' TODO Would be convenient to handle multiple samples, as sample1/xxx; in this case, 
#' should prepend the sample name to the barcodes.
#' 
#' issue: when shardifying, good to keep info about what to merge. this reduces the work plenty!
#' could keep a list of which shards belong for the next step
#' 
#' 
#' @param outputFiles Files that we expect to exist
#' @param overwrite Files that we expect to exist
#' @return TRUE if ok to proceed
bascet_check_overwrite_output <- function(
  outputFiles, 
  overwrite
){
  if(all(file.exists(outputFiles)) & !overwrite){
    print("All files to write already exist; skipping. To change this behaviour, set overwrite=TRUE")
    FALSE
  } else {
    TRUE
  }
}

  


################################################################################
################ Bascet command line tool: debarcoding raw fastq ###############
################################################################################



###############################################
#' Detect metadata for raw input FASTQ files
#' 
#' _R1 -- common from illumina sequencer
#' SRR****_1.fastq.gz -- typical from SRA
#' 
#' TODO Would be convenient to handle multiple samples, as sample1/xxx; in this case, 
#' should prepend the sample name to the barcodes.
#' 
#' issue: when shardifying, good to keep info about what to merge. this reduces the work plenty!
#' could keep a list of which shards belong for the next step
#' 
#' 
#' @param rawRoot Path to folder with FASTQ files
#' @export
#' @return A data frame with metadata for the raw input files
DetectRawFileMeta <- function(rawRoot, verbose=FALSE){
  #rawRoot <- "/husky/fromsequencer/241206_novaseq_wgs3/raw"
  allfiles <- list.files(rawRoot)
  
  if(verbose){
    print("Found the following files in the directory")
    print(allfiles)
  }
  
  allfiles <- allfiles[
    stringr::str_ends(allfiles,".fq.gz") | 
      stringr::str_ends(allfiles, ".fastq.gz") | 
      stringr::str_ends(allfiles, ".fq") | 
      stringr::str_ends(allfiles,".fastq")]
  
  if(verbose){
    print("Found the following raw read-like files")
    print(allfiles)
  }
  
  r1_files <- allfiles[stringr::str_detect(allfiles,stringr::fixed("_R1"))]
  
  if(verbose){
    print("Found the following R1-like files")
    print(r1_files)
  }
  
  ## TODO for SRA, replace with _1.fastq.gz
  r2_corresponding <- stringr::str_replace(r1_files,stringr::fixed("R1"),"R2")
  
  if(!all(r2_corresponding %in% allfiles)){
    stop("Not all R1 files appear to have a corresponding R2 file")
  }
  
  ### Gather first set of metadata
  meta <- data.frame(
    r1=r1_files,
    r2=r2_corresponding,
    dir=file.path(rawRoot) #this standarizes the trailing /
  )
  meta$prefix <- ""
  
  #Guess prefixes
  meta$possible_prefix <- stringr::str_split_i(meta$r1, "_S[0123456789]",1)
  unique_prefix <- unique(meta$possible_prefix)
  if(verbose){
    print("Detected possible prefixes")
    print(unique_prefix)
  }

  #If there is more than one prefix, then we have to add them. otherwise just keep it simple
  if(length(unique_prefix)>1){
    #Sanitize prefixes. Some characters will break BAM tags etc
    meta$prefix <- stringr::str_remove_all(meta$prefix, " ")
    meta$prefix <- stringr::str_remove_all(meta$prefix, "/")
    meta$prefix <- stringr::str_remove_all(meta$prefix, "\"")
    
  } else {
    if(verbose){
      print("Detected one one possible prefix, so not adding it")
    }
  }

  #Return metadata
  meta[,c("prefix","r1","r2","dir")]
}

#rawmeta <- DetectRawFileMeta("/husky/fromsequencer/241206_novaseq_wgs3/raw")


###############################################
#' Generate BAM with barcodes from input raw FASTQ
#' 
#' @inheritParams template_BascetFunction
#' @param rawmeta Metadata for the raw FASTQ input files. See DetectRawFileMeta
#' @param outputName Name output files: Debarcoded reads
#' @param outputNameIncomplete Name of output files: Reads that could not be parsed
#' @param chemistry The type of data to be parsed
#' @param overwrite description
#' @export
BascetGetRaw <- function(
    bascetRoot, 
    rawmeta, 
    outputName="debarcoded", 
    outputNameIncomplete="incomplete_reads", 
    chemistry="atrandi_wgs",  #or atrandi_rnaseq; any way to get list from software?
    overwrite=FALSE,
    runner, 
    bascet_instance=bascet_instance.default
){

  #One output per input pair of reads  
  num_shards <- nrow(rawmeta)
  if(num_shards==0){
    stop("No input files")
  }
  outputFilesComplete <- make_output_shard_names(bascetRoot, outputName, "tirp.gz", num_shards)
  outputFilesIncomplete <- make_output_shard_names(bascetRoot, outputNameIncomplete, "tirp.gz", num_shards)
  
  #Check if libnames should be added
  add_libnames <- any(rawmeta$prefix!="")

  if(bascet_check_overwrite_output(outputFilesComplete)) {
    RunJob(
      runner = runner, 
      jobname = "bascet_getraw",
      cmd = c(
        shellscript_set_tempdir(bascet_instance),
        shellscript_make_bash_array("files_r1",file.path(rawmeta$dir, rawmeta$r1)),
        shellscript_make_bash_array("files_r2",file.path(rawmeta$dir, rawmeta$r2)),
        shellscript_make_bash_array("libnames",rawmeta$prefix),
        shellscript_make_bash_array("files_out",outputFilesComplete),
        shellscript_make_bash_array("files_out_incomplete",outputFilesIncomplete),
        paste(
          bascet_instance@prepend_cmd,
          bascet_instance@bin, 
          "getraw",
          "-t $BASCET_TEMPDIR",
          "--chemistry",chemistry,  
          "--r1 ${files_r1[$TASK_ID]}",  
          "--r2 ${files_r2[$TASK_ID]}",
          if(add_libnames) "--libname ${libnames[$TASK_ID]}",
          "--out-complete   ${files_out[$TASK_ID]}",                 #Each job produces a single output
          "--out-incomplete ${files_out_incomplete[$TASK_ID]}"                 #Each job produces a single output
        )
      ),
      arraysize = nrow(rawmeta)
    )    
  } else {
    new_no_job()
  }
}





################################################################################
################ Bascet command line tools: WGS ################################
################################################################################




###############################################
#' Take debarcoded reads and split them into suitable numbers of shards.
#' 
#' The reads from one cell is guaranteed to only be present in a single shard.
#' This makes parallel processing simple as each shard can be processed on
#' a separate computer. Using more shards means that more computers can be used.
#' 
#' If you perform all the calculations on a single computer, having more
#' than one shard will not result in a speedup. This option is only relevant
#' when using a cluster of compute nodes.
#' 
#' 
#' 
#' TODO if we have multiple input samples, is there a way to group them?
#' otherwise we will be reading more input files than needed. that said,
#' if we got an index, so if list of cells specified, it is possible to quickly figure out
#' out if a file is needed at all for a merge
#' 
#' TODO Figuring out if a file is needed can be done at "planning" (Zorn) stage
#' 
#' TODO seems faster to have a single merger that writes multiple output files if
#' cell list is not provided. if the overhead is accepted then read all input files and
#' discard cells on the fly
#' 
#' @param inputName Name of input file: Debarcoded reads
#' @param outputName Name of the output file: Properly sharded debarcoded reads
#' 
#' @export
BascetShardify <- function(
    bascetRoot, 
    inputName="debarcoded", 
    includeCells=NULL, ############# TODO: get rid of this parameter; only support direct merging
    num_output_shards=1,
    outputName="filtered", 
    overwrite=FALSE,
    runner, 
    bascet_instance=bascet_instance.default
){

  input_shards <- detect_shards_for_file(bascetRoot, inputName)
  if(length(input_shards)==0){
    stop("Found no input files")
  }

  #Include all cells if nothing else provided
  if(is.null(includeCells)){
    includeCells <- BascetCellNames(bascetRoot, inputName)
    includeCells <- unique(includeCells$cell) #when shardifying, we expect cells to appear more than once  -- could warn for other commands!
    print(paste("Including all the",length(includeCells), "cells"))
  }

  #Figure out which cell goes into which shard
  list_cell_for_shard <- shellscript_split_arr_into_list_randomly(includeCells, num_output_shards)
  
  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, input_shards)
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "tirp.gz", num_output_shards)

  if(bascet_check_overwrite_output(outputFiles, overwrite)) {
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "bascet_shardify",
      cmd = c(
        shellscript_set_tempdir(bascet_instance),
        shellscript_make_files_expander("CELLFILE", list_cell_for_shard),
        shellscript_make_bash_array("files_out", outputFiles),
        paste(
          bascet_instance@prepend_cmd,
          bascet_instance@bin,
          "shardify", 
          "-t $BASCET_TEMPDIR",
          "-i",shellscript_make_commalist(inputFiles), #Need to give all input files for each job
          "-o ${files_out[$TASK_ID]}",                 #Each job produces a single output
          "--cells $CELLFILE"                          #Each job takes its own list of cells
        ),  
        "rm $CELLFILE"
      ), 
      arraysize = num_output_shards
    )
  } else {
    new_no_job()
  }
}




###############################################
#' Transform data
#' 
#' This command enables
#' * subsetting to a list of cells
#' * converting between file formats
#' * merging shards
#' * dividing shards
#' 
#' @param num_divide description
#' @param num_merge description
#' @param out_format description
#' 
#' @export
BascetMapTransform <- function(
    bascetRoot, 
    inputName, 
    outputName,
    num_divide=1,
    num_merge=1,
    out_format="tirp.gz", ### not really!
    includeCells=NULL,
    overwrite=FALSE,
    runner, 
    bascet_instance=bascet_instance.default
){
  
  
  is.integer.overequal1 <- function(x){
    return(x>=1 & x==as.integer(x))
  }
  
  #Verify input
  if(!(num_divide==1 || num_merge==1)){
    stop("Either divide or merge must be set to 1")
  }
  if(!is.integer.overequal1(num_divide) || !is.integer.overequal1(num_merge)){
    stop("Number of divisions and merges must be an integer >= 1")
  }
  
  #TODO check outformat

  #Figure out input and output file names  
  inputFiles <- file.path(bascetRoot, detect_shards_for_file(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }
  
  outputFiles <- make_output_shard_names(bascetRoot, outputName, out_format, num_shards)

  #If cell list is provided, produce a file for input (not all transform calls can handle this, so optional)
  produce_cell_list <- !is.null(includeCells)
  if(produce_cell_list) {
    #Figure out which cell goes into which shard   TODO, using the same cell list for all shards
    list_cell_for_shard <- list()
    for(i in 1:length(inputFiles)){
      list_cell_for_shard[[i]] <- includeCells
    }
  }
  
  if(bascet_check_overwrite_output(outputFiles, overwrite)) {
    #Run the job
    RunJob(
      runner = runner, 
      jobname = paste0("bascet_transform"),
      cmd = c(
        shellscript_set_tempdir(bascet_instance),
        if(produce_cell_list) shellscript_make_files_expander("CELLFILE", list_cell_for_shard),
        shellscript_make_bash_array("files_in",inputFiles),
        shellscript_make_bash_array("files_out",outputFiles),
        paste(
          bascet_instance@prepend_cmd,
          bascet_instance@bin, 
          "transform",
          if(produce_cell_list) "--cells $CELLFILE",
          #"-t $BASCET_TEMPDIR",  ##not supported
          "-i ${files_in[$TASK_ID]}",
          "-o ${files_out[$TASK_ID]}"
        )
      ),
      arraysize = num_shards
    )  
  } else {
    new_no_job()
  }
}






################################################################################
################ Quality control ###############################################
################################################################################



###############################################
#' 
#' Produce a kneeplot
#' 
#' @param adata description
#' @param max_species description
#' 
#' @return A ggplot object
#' @export
KneeplotPerSpecies <- function(
    adata, 
    max_species=NULL
) {
  strain_cnt <- adata@assays[[DefaultAssay(adata)]]$counts
  
  if(!is.null(max_species)){
    strain_cnt <- strain_cnt[order(rowSums(strain_cnt), decreasing = TRUE),]
    strain_cnt <- strain_cnt[1:min(nrow(strain_cnt), max_species),]
  }
  
  allknee <- list()
  for(i in 1:nrow(strain_cnt)){
    onedf <- data.frame(
      strain = rownames(strain_cnt)[i],
      cnt = strain_cnt[i,]
    )
    onedf <- onedf[order(onedf$cnt, decreasing = TRUE),]
    onedf$index <- 1:nrow(onedf)
    allknee[[paste("s",i)]] <- onedf
  }
  allknee <- do.call(rbind, allknee)
  
  ggplot(allknee, aes(index, cnt, color=strain)) + geom_line() +
    scale_x_log10() +
    scale_y_log10() +
    xlab("Log10 Cell index") +
    ylab("Log10 Total read count") +
    theme_bw()
  # +
  #  theme(legend.position = "none")
  
}




###############################################
#' 
#' Produce a matrix of Barnyard plots, i.e., counts for one species vs another, 
#' for all combinations of species.
#' 
#' 
#' 
#' 
#' @param adata A Seurat object with the DefaultAssay set to species count
#' @return a ggarranged set of ggplots
#' @export
BarnyardPlotMatrix <- function(
    adata
){
  cnt <- adata@assays[[DefaultAssay(adata)]]$counts
  cnt <- cnt[rowSums(cnt)>0,]  ##Only consider species we have
  list_species <- rownames(cnt)#[1:5] ######## Just do a few!
  
  all_plots <- list()
  for(i in seq_along(list_species)){
    for(j in seq_along(list_species)){
      #print(paste(i,j))
      if(i>=j) {
        all_plots[[paste(i,j)]] <- ggplot()
      } else {
        df <- data.frame(
          x=cnt[i,],
          y=cnt[j,]
        )
        p <- ggplot(df, aes(x+1,y+1)) +
          geom_point() + 
          scale_x_log10() + 
          scale_y_log10() + 
          theme_bw()+
          xlab(paste("Pseudocount", list_species[i]))+
          ylab(paste("Pseudocount", list_species[j]))
        
        all_plots[[paste(i,j)]] <- p        
      }
    }
  }
  egg::ggarrange(plots = all_plots, nrow = length(list_species))  
}




