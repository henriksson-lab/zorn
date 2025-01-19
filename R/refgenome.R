 


################################################################################
################## Alignment ###################################################
################################################################################





###############################################
#' Align from FASTQ, generate sorted and indexed BAM file
#' 
#' TODO Should be managed by Bascet mapshard system, with automatic input conversion. unaligned file should be made temp and removed
#' 
#' @return TODO
#' @export
BascetAlignToReference <- function(
    bascetRoot, 
    useReference,
    numLocalThreads=1,
    inputName="asfq", ######### should be able to take filtered and pipe to bwa if needed  "filtered"
    outputNameBAMunsorted="unsorted_aligned", 
    outputNameBAMsorted="aligned", 
    runner, 
    bascet_instance=bascet_instance.default){
  
  #Figure out input and output file names  
  input_shards <- detect_shards_for_file(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards) #### TODO 666 should really need to add this?
  
  outputFilesBAMunsorted <- make_output_shard_names(bascetRoot, outputNameBAMunsorted, "bam", num_shards)
  outputFilesBAMsorted   <- make_output_shard_names(bascetRoot, outputNameBAMsorted,   "bam", num_shards)
  
  #### Align with BWA
  cmd <- paste(
    "bwa mem",
    useReference,
    inputFiles[1],
    
#    file.path(bascetRoot,"for_bwa.1.fq.gz"),
    "-t ",numLocalThreads,
    "> ",
    outputFilesBAMunsorted[1]
  )
  system(cmd)
  
  #### Sort the aligned reads
  cmd <- paste("samtools sort",
               outputFilesBAMunsorted[1],
               "-@ ",numLocalThreads,
               "-o ",
               outputFilesBAMsorted[1]
  )
  system(cmd)
  
  #### Index the file
  cmd <- paste("samtools index",
               outputFilesBAMsorted[1],
               "-@ ",numLocalThreads
  )
  system(cmd)
  
  ### TODO this should be a proper job!
  
  666
}









###############################################
#' Take aligned BAM file and produce Fragments.tsv.gz, compatible with Signac ATAC-seq style analysis
#' 
#' @return TODO
#' @export
BascetBam2Fragments <- function(
    bascetRoot, 
    inputName="aligned",
    outputName="fragments", 
    runner, 
    bascet_instance=bascet_instance.default){
  
  input_shards <- detect_shards_for_file(bascetRoot, inputName)
  if(length(input_shards)==0){
    stop("Found no input files")
  }
  num_shards <- length(input_shards)
  inputFiles <- file.path(bascetRoot, input_shards)
  
  #One output per input pair of reads  
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "tsv.gz", num_shards)  
  

  RunJob(
    runner = runner, 
    jobname = "bascet_bam2fragments",
    cmd = c(
      shellscript_set_tempdir(bascet_instance),
      shellscript_make_bash_array("files_in", inputFiles),
      shellscript_make_bash_array("files_out",outputFiles),
      paste(
        bascet_instance@bin, 
        "bam2fragments",
        "-t $BASCET_TEMPDIR",
        "-i ${files_in[$TASK_ID]}",  
        "-o ${files_out[$TASK_ID]}"
      )
    ),
    arraysize = nrow(rawmeta)
  )
}









###############################################
#' Using Tabix, get list of sequences in a fragment file
#' 
#' @return TODO
#' @export
TabixGetFragmentsSeqs <- function(fragpath) {
  system(paste("tabix -l ", fragpath), intern = TRUE)
}


###############################################
#' From a fragments file, get a chromatin assay for Signac
#' 
#' @return TODO
#' @export
FragmentsToSignac <- function(fragpath) {
  #### Index if needed; this should already have been done but added here just in case
  fragpath_index <- paste(fragpath,".tbi",sep="")
  if(!file.exists(fragpath_index)){
    print("Indexing fragment file")
    system(paste("tabix -p bed ",fragpath))
  }
  
  #### Count fragments / see which cells are present. Possible to add a cutoff here if we wish
  counts <- CountFragments(fragments = fragpath)  
  colnames(counts)[1] <- "cb"
  counts <- counts[order(counts$reads_count, decreasing = TRUE),]
  
  #knee plot analysis
  if(FALSE){
    kneedf <- data.frame(cnt=sort(counts$reads_count, decreasing = TRUE))
    kneedf$index <- 1:nrow(kneedf)
    ggplot(kneedf, aes(index, cnt)) + 
      geom_line() + 
      scale_x_log10() + 
      scale_y_log10()
  }
  
  keep_cells <- counts$cb #[log10(counts$reads_count)>0]
  
  #### Create a dummy assay
  stupidmat <- matrix(0, nrow = 2, ncol=length(keep_cells))
  rownames(stupidmat) <- c("chr1:1-100", "chr1:200-300")
  colnames(stupidmat) <- keep_cells
  
  chrom_assay <- CreateChromatinAssay(
    counts = as.sparse(stupidmat),
    sep = c(":", "-"),
    fragments = fragpath,
    min.cells = 0,
    min.features = 0
  )  
  chrom_assay
}



###############################################
#' From a Signac chromatin assay with fragments, for each cell, count how many reads per chromosome
#' 
#' @return TODO
#' @export
FragmentCountsPerChrom <- function(chrom_assay){
  
  fr <- Fragments(chrom_assay)[[1]]
  #fr@cells ## also possible!
  
  all_seqid <- TabixGetFragmentsSeqs(fr@path)
  
  grange <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqid=all_seqid,
    name=all_seqid,
    start=1,
    end=536870912  #for grange this must be <2**31. for scanTabix, this must be <=536870912
  ))
  
  # quantify counts over all chrom
  reg_counts <- FeatureMatrix(
    fragments = Fragments(chrom_assay),
    features = grange,
    cells = colnames(fr@cells)
  )  
  
  reg_counts
}




###############################################
#' From a Signac chromatin assay with fragments, for each cell, count how many reads per chromosome.
#' This function directly returns an assay that can be added to a Seurat multimodal object
#' 
#' @return TODO
#' @export
FragmentCountsPerChromAssay <- function(
    bascetRoot,
    inputName="fragments.1.tsv.gz" ### TODO, detect, merge all of them
) {
  fragpath <- file.path(bascetRoot,"fragments.1.tsv.gz")
  chrom_assay <- FragmentsToSignac(fragpath)
  chrom_assay <- FragmentCountsPerChrom(chrom_assay)
  
  #Subset; needed?
  #include_cells <- sort(intersect(colnames(adata), colnames(chrom_assay)))
  #chrom_assay <- chrom_assay[,include_cells]
  #adata <- adata[,include_cells]
  
  chrom_assay_hack <- CreateAssay5Object(counts = chrom_assay)
  rownames(chrom_assay_hack) <- stringr::str_split_i(rownames(chrom_assay_hack),"-1-",1)  #evil hack, fragile

  #Figure out which chromosome dominates  
  #chrom_assay_hack@meta.data$max_chrom <- rownames(chrom_assay_hack$counts)[apply(chrom_assay_hack$counts, 2, which.max)]
  
  chrom_assay_hack
  #obj <- CreateAssay5Object(counts = chrom_assay_hack)
  #obj$max_chrom <- max_chrom
  #obj
}




###############################################
#' Produce a count matrix on strain level. too specific?
#' @return TODO
#' @export
ChromToSpeciesCount <- function(adata, map_seq2strain){
  mat_cnt <- adata@assays[[DefaultAssay(adata)]]$counts
  
  unique_strains <- unique(map_seq2strain$strain)
  strain_cnt <- matrix(NA, ncol=ncol(mat_cnt), nrow=length(unique_strains))
  rownames(strain_cnt) <- unique_strains
  colnames(strain_cnt) <- colnames(mat_cnt)
  for(i in 1:nrow(strain_cnt)){
    cur_cols <- map_seq2strain$id[map_seq2strain$strain %in% rownames(strain_cnt)[i]]
    print(cur_cols)
    strain_cnt[i,] <- colSums(mat_cnt[rownames(mat_cnt) %in% cur_cols,,drop=FALSE])
  }
  
  CreateAssay5Object(counts = strain_cnt)
}


