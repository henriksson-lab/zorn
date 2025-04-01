 


################################################################################
################## Alignment ###################################################
################################################################################




###############################################
#' Generate a bigwig out of all reads in a sorted BAM. Note that the caller
#' is responsible for sorting the BAM first
#' 
#' This function is mainly for QC purposes. It uses bamCoverage from deepTools
#' apt install python3-deeptools
#' 
#' @param outputName Name of output file (BIGWIG-file)
#' 
BascetAlignmentToBigwig <- function(  
    bascetRoot, 
    inputName="aligned", 
    outputName="pileup",
    runner, 
    bascet_instance=bascet_instance.default){
  
  #Figure out input and output file names
  input_shards <- detect_shards_for_file(bascetRoot, inputName) 
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards) 
  
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "bigwig", num_shards)
  
  
  ### Build command
  cmd <- c(
    shellscript_set_tempdir(bascet_instance),
    shellscript_make_bash_array("files_in", inputFiles),
    shellscript_make_bash_array("files_out", outputFiles),
    
    ### For sorting
    paste(
      bascet_instance@prepend_cmd,
      "bamCoverage",
      "-b ${files_in[$TASK_ID]}",     #Each job takes a single output
      "-o ${files_out[$TASK_ID]}"     #Each job produces a single output
    )
  )
  
  
  
  #Produce the script and run the job
  RunJob(
    runner = runner, 
    jobname = "bascet_tobw",
    cmd = cmd,
    arraysize = num_shards
  )
}










###############################################
#' Filter an alignment (BAM-file).
#' 
#' This is typically used to either remove host DNA, or keep reads mapping to a known reference
#' 
#' @param outputName Name of output file (BAM-file)
#' @param numLocalThreads Number of threads to use for each runner
#' @param keep_mapped Keep the mapped reads (TRUE) or unmapped (FALSE)
#' 
BascetFilterAlignment <- function(  
    bascetRoot, 
    numLocalThreads=1,
    inputName, 
    outputName,
    keep_mapped=FALSE,
    runner, 
    bascet_instance=bascet_instance.default){
  
  #Figure out input and output file names
  input_shards <- detect_shards_for_file(bascetRoot, inputName) 
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards) 
  
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "bam", num_shards)

  ### Set flags for samtools
  samtools_flags <- ""
  if(keep_mapped){
    samtools_flags <- paste(samtools_flags, "-F4")
  } else {
    samtools_flags <- paste(samtools_flags, "-f4")
  }
  
  ### Build command
  cmd <- c(
    shellscript_set_tempdir(bascet_instance),
    shellscript_make_bash_array("files_in", inputFiles),
    shellscript_make_bash_array("files_out", outputFiles),
    
    ### For sorting
    paste(
      bascet_instance@prepend_cmd,
      "samtools view",
      samtools_flags, 
      "-@ ",numLocalThreads,          #Number of threads to use
      "${files_in[$TASK_ID]}",        #Each job takes a single output
      "-o ${files_out[$TASK_ID]}"     #Each job produces a single output
    )
  )
    
  
  
  #Produce the script and run the job
  RunJob(
    runner = runner, 
    jobname = "bascet_filteraln",
    cmd = cmd,
    arraysize = num_shards
  )
}






## override template for docs


###############################################
#' Align from FASTQ, generate sorted and indexed BAM file
#' 
#' TODO Should be managed by Bascet mapshard system, with automatic input conversion. unaligned file should be made temp and removed
#' 
#' @param useReference Name of the BWA reference to use
#' @param outputNameBAMunsorted Name of unsorted BAMs
#' @param outputNameBAMsorted Name of sorted BAMs (if generated)
#' @param numLocalThreads Number of threads to use for each runner
#' @param do_sort Whether to sort the output or not
#' 
#' @export
BascetAlignToReference <- function(  
    bascetRoot, 
    useReference,
    numLocalThreads=1,
    inputName="asfq", ######### should be able to take filtered and pipe to bwa if needed  "filtered"
    outputNameBAMunsorted="unsorted_aligned", 
    outputNameBAMsorted="aligned", 
    do_sort=TRUE,
    runner, 
    bascet_instance=bascet_instance.default){
  
  #Figure out input and output file names  
  input_shards <- detect_shards_for_file(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards) 
  
  outputFilesBAMunsorted <- make_output_shard_names(bascetRoot, outputNameBAMunsorted, "bam", num_shards)
  outputFilesBAMsorted   <- make_output_shard_names(bascetRoot, outputNameBAMsorted,   "bam", num_shards)
  
  ### Build command: basic alignment
  cmd <- c(
    shellscript_set_tempdir(bascet_instance),
    shellscript_make_bash_array("files_in", inputFiles),
    shellscript_make_bash_array("files_out_unsorted", outputFilesBAMunsorted),
    shellscript_make_bash_array("files_out_sorted", outputFilesBAMsorted),
    
    ### For alignment
    paste(
      bascet_instance@prepend_cmd,
      "bwa mem", 
      useReference,
      "${files_in[$TASK_ID]}",                #Align one input
      "-t", numLocalThreads,
      ">",
      "${files_out_unsorted[$TASK_ID]}"       #Each input means one output
    )
  )
  
  
  ### Build command: sorting and indexing
  if(do_sort){
    cmd <- c(
      cmd,
      
      ### For sorting
      paste(
        bascet_instance@prepend_cmd,
        "samtools sort", 
        "-@ ",numLocalThreads,                 #Number of threads to use
        #        "-T $BASCET_TEMPDIR",                          #temporary FILE. we only got directory... TODO
        "${files_out_unsorted[$TASK_ID]}",    #Each job produces a single output
        "-o ${files_out_sorted[$TASK_ID]}"     #Each job produces a single output
      ),
      
      ### For indexing
      paste(
        bascet_instance@prepend_cmd,
        "samtools index", 
        "-@ ",numLocalThreads,                 #Number of threads to use
        "${files_out_sorted[$TASK_ID]}"        #Each job produces a single output
      )      
    )
  }
  
  
  #Produce the script and run the job
  RunJob(
    runner = runner, 
    jobname = "bascet_aln",
    cmd = cmd,
    arraysize = num_shards
  )
}









###############################################
#' Take aligned BAM file and produce Fragments.tsv.gz, compatible with Signac ATAC-seq style analysis
#' 
#' @return A job, producing a type of Fragments.tsv.gz
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
        bascet_instance@prepend_cmd,
        bascet_instance@bin, 
        "bam2fragments",
        "-t $BASCET_TEMPDIR",
        "-i ${files_in[$TASK_ID]}",  
        "-o ${files_out[$TASK_ID]}"
      )
    ),
    arraysize = length(inputFiles)
  )
}





###############################################
#' From aligned BAM file, compute counts per chromosome
#' 
#' @return A job, executing the counting 
#' @export
BascetCountChrom <- function(
    bascetRoot, 
    inputName="aligned",
    outputName="chromcount", 
    runner, 
    bascet_instance=bascet_instance.default){
  
  input_shards <- detect_shards_for_file(bascetRoot, inputName)
  if(length(input_shards)==0){
    stop("Found no input files")
  }
  num_shards <- length(input_shards)
  inputFiles <- file.path(bascetRoot, input_shards)
  
  #One output per input alignment
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "hd5", num_shards)  
  
  
  RunJob(
    runner = runner, 
    jobname = "bascet_countchrom",
    cmd = c(
      shellscript_set_tempdir(bascet_instance),
      shellscript_make_bash_array("files_in", inputFiles),
      shellscript_make_bash_array("files_out",outputFiles),
      paste(
        bascet_instance@prepend_cmd,
        bascet_instance@bin, 
        "countchrom",
        "-t $BASCET_TEMPDIR",
        "-i ${files_in[$TASK_ID]}",  
        "-o ${files_out[$TASK_ID]}"
      )
    ),
    arraysize = length(inputFiles)
  )
}





###############################################
#' Using Tabix, get list of sequences in a fragment file
#' 
#' TODO too raw? wrap in bascet?
#' 
#' @return The raw result of tabix
#' @export
TabixGetFragmentsSeqs <- function(
    fragpath,
    bascet_instance=bascet_instance.default
) {
  system(paste(bascet_instance@prepend_cmd,"tabix -l ", fragpath), intern = TRUE)
}


###############################################
#' From a fragments file, get a chromatin assay for Signac  .......................... signac is dirt slow at counting! only 1 cpu used. take over the whole process
#' 
#' Case: tabix from command line uses 100% cpu only. likely not designed for large queries
#' 
#' 
#' @return a ChromatinAssay
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
#' @return a FeatureMatrix 
#' @export
FragmentCountsPerChrom <- function(
    chrom_assay,
    bascet_instance=bascet_instance.default
){
  
  #Figure out where the fragment file is
  fr <- Fragments(chrom_assay)[[1]]
  #fr@cells ## also possible!
  
  #Get the name of chromosomes
  all_seqid <- TabixGetFragmentsSeqs(fr@path, bascet_instance)
  grange <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqid=all_seqid,
    name=all_seqid,
    start=1,
    end=536870912  #for grange this must be <2**31. for scanTabix, this must be <=536870912
  ))
  
  # Quantify counts over each chromosome
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
#' @return A seurat object holding the counts
#' @export
FragmentCountsPerChromAssay <- function(
    bascetRoot,
    inputName="fragments.1.tsv.gz" ### TODO, detect, merge all of them
) {
  fragpath <- file.path(bascetRoot,inputName)
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
#' Produce a count matrix on strain level
#' 
#' TODO too specific?
#' 
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






################################################################################
############ RNA-seq style feature counting from fragments.tsv #################
################################################################################


###############################################
#' Obtain a feature matrix (as seurat object) given an seurat object having Fragments associated
#' 
#' @param adata Seurat object, with FeatureMatrix
#' @param grange_gene A grange object telling where to count
#' @return A Seurat AssayObject
#' @export
CountGrangeFeatures <- function(adata, grange_gene){
  gene_counts <- FeatureMatrix(
    fragments = Fragments(adata),
    features = grange_gene,
    cells = colnames(adata)
  )
  rownames(gene_counts) <- grange_gene$Name ## paste(grange_gene$gene_biotype,grange_gene$Name) #Name in an easy manner
  
  # create a new assay using the MACS2 peak set 
  CreateAssayObject(
    counts = gene_counts[!duplicated(rownames(gene_counts)),],
    min.cells = 0,  #filter later!!!
    min.features = 0
  )
}


