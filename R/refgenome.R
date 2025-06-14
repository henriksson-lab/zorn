 

################################################################################
################## Alignment ###################################################
################################################################################


###############################################
#' Figure out if a BAM-file is a paired alignment or not
#' 
#' @param fname Name of BAM-file
#' 
#' @export
is_bam_paired_alignment <- function(
    fname,
    bascet_instance=GetDefaultBascetInstance()
){

  #example
  #@PG     ID:bwa  PN:bwa  VN:0.7.18-r1243-dirty   CL:bwa mem /home/m/mahogny/mystore/atrandi/bwa_ref/combo_hx/all.fa.gz /pfs/proj/nobackup/fs/projnb10/hpc2nstor2024-027/atrandi/v2_wgs_saliva1/asfq.1.R1.fq.gz /pfs/proj/nobackup/fs/projnb10/hpc2nstor2024-027/atrandi/v2_wgs_saliva1/asfq.1.R2.fq.gz -t 10
  
  #Note: can also see this from the BAM header, but a real mess to parse out!
  # samtools view unsorted_aligned.1.bam  -H
  # @PG	ID:bwa	PN:bwa	VN:0.7.17-r1198-dirty	CL:bwa mem /husky/fromsequencer/241210_joram_rnaseq/ref/all.fa /husky/henriksson/atrandi/v2_rnaseq5//asfq.1.R1.fq.gz /husky/henriksson/atrandi/v2_rnaseq5/asfq.1.R2.fq.gz -t 20
  
  cmd <- paste(
    bascet_instance@prepend_cmd, 
    "samtools view",
    fname,
    "-h | head -n 10000 | ", #100 was not enough to be reliable
    bascet_instance@prepend_cmd," samtools view -c -f 1"
  )
  
  ret <- system(
    cmd, 
    intern = TRUE, 
    ignore.stderr=TRUE,
  )
  
  print("Ignore samtools view errors here; this should be rewritten to be handled by bascet instead")
  as.integer(ret)>0
}




###############################################
#' Index a genome using BWA such that it can be used for alignment
#' 
#' @param genomeFile Name of FASTA file
#' 
#' @export
BascetIndexGenomeBWA <- function(  
    bascetRoot, 
    genomeFile, 
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascet_instance=GetDefaultBascetInstance()
){
  
  if(!file.exists(genomeFile)){
    stop("Could not find genome FASTA file")
  }
  
  if(bascet_check_overwrite_output(outputFiles, overwrite)) {
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Z_bwa_index",
      cmd = c(
        ### For sorting
        paste(
          bascet_instance@prepend_cmd,
          "bwa index",
          genomeFile
        )
      ),
      arraysize = 1
    )
  } else {
    new_no_job()
  }
}



###############################################
#' Index a genome using STAR such that it can be used for alignment
#' 
#' @param genomeFile Name of FASTA file
#' @param gtfFile description
#' @param outDir A directory in which to store the output genome. This directory will be created
#' 
#' @export
BascetIndexGenomeSTAR <- function(  
    bascetRoot, 
    fastaFile, 
    gtfFile,
    outDir,
    numLocalThreads=10,
    runner=GetDefaultBascetRunner(), 
    bascet_instance=GetDefaultBascetInstance()
){
  
  if(!file.exists(fastaFile)){
    stop("Could not find genome FASTA file")
  }
  if(!file.exists(gtfFile)){
    stop("Could not find genome GTF file")
  }
  if(file.exists(outDir)){
    return(new_no_job())
  }

    
  #Produce the script and run the job.
  #Note: overwrite not supported. It is too dangerous to implement
  RunJob(
    runner = runner, 
    jobname = "Z_star_index",
    cmd = c(
      paste(
        "mdkir -p", outDir
      ),
      paste(
        bascet_instance@prepend_cmd,
        "STAR --runThreadN", as.character(numThreads),
        "--runMode genomeGenerate",
        "--genomeDir",outDir,
        "--genomeFastaFiles",fastaFile,
        "--sjdbGTFfile", gtfFile
      )
    ),
    arraysize = 1
  )
}






###############################################
#' Check if a file is a FASTQ file
#' 
#' @param fname Path to file
#' @return TRUE if the file is some type of FASTQ
is_fastq <- function(fname) {
  stringr::str_ends(fname, stringr::fixed("fq.gz")) ||  
    stringr::str_ends(fname, stringr::fixed("fastq.gz"))
}

###############################################
#' Check if a file is a paired FASTQ file
#' 
#' Panics if the file is not a FASTQ at all
#' 
#' @param fname Path to file
#' @return TRUE if the file is a paired FASTQ
is_paired_fastq <- function(fname) {
  if(is_fastq(fname)){
    stringr::str_ends(fname, stringr::fixed("R1.fq.gz")) ||  
      stringr::str_ends(fname, stringr::fixed("R1.fastq.gz")) ||
      stringr::str_ends(fname, stringr::fixed("R2.fq.gz")) ||  
      stringr::str_ends(fname, stringr::fixed("R2.fastq.gz"))
  } else {
    stop("The file is not any type of FASTQ")
  }
}



###############################################
#' Get corresponding R2 file. Assumes that the input file is R1
#' 
#' @param fname Path to file
#' @return Path to R2 file
get_fastq_R2_from_R1 <- function(fname) {

  dir <- dirname(fname)
  bname <- basename(fname)

  bname <- stringr::str_replace(bname, stringr::fixed("R1.fq.gz"), "R2.fq.gz")  
  bname <- stringr::str_replace(bname, stringr::fixed("R1.fastq.gz"), "R2.fastq.gz")  

  file.path(dir, bname)
}





###############################################
#' Generate a bigwig out of all reads in a sorted BAM. Note that the caller
#' is responsible for sorting the BAM first
#' 
#' This function is mainly for QC purposes. It uses bamCoverage from deepTools
#' apt install python3-deeptools
#' 
#' @param outputName Name of output file (BIGWIG-file)
#' 
#' @export
BascetAlignmentToBigwig <- function(  
    bascetRoot, 
    inputName="aligned", 
    outputName="pileup",
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascet_instance=GetDefaultBascetInstance()
){
  
  #Figure out input and output file names
  input_shards <- detect_shards_for_file(bascetRoot, inputName) 
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards) 
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "bigwig", num_shards)
  
  if(bascet_check_overwrite_output(outputFiles, overwrite)) {
    
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
      jobname = "Z_tobw",
      cmd = cmd,
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
}










###############################################
#' Filter an alignment (BAM-file).
#' 
#' This is typically used to either remove host DNA, or keep reads mapping to a known reference
#' 
#' If the BAM-file has paired reads then BOTH reads need to be mapped (flag 0x2); otherwise (flag 0x4)
#' 
#' 
#' @param outputName Name of output file (BAM-file)
#' @param numLocalThreads Number of threads to use for each runner
#' @param keep_mapped Keep the mapped reads (TRUE) or unmapped (FALSE)
#' 
#' @export
BascetFilterAlignment <- function(  
    bascetRoot, 
    numLocalThreads=1,
    inputName, 
    outputName,
    keep_mapped=FALSE,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascet_instance=GetDefaultBascetInstance()
){
  
  
  #Figure out input and output file names
  input_shards <- detect_shards_for_file(bascetRoot, inputName) 
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards) 
  
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "bam", num_shards)

  if(bascet_check_overwrite_output(outputFiles, overwrite)) {
    
    ### Set flags for samtools. Depends on if paired alignment or not
    samtools_flags <- ""
    is_paired_al <- is_bam_paired_alignment(inputFiles[1], bascet_instance)
    print(paste("Detect paired alignment: ",is_paired_al))
    
    if(is_paired_al) {
      if(keep_mapped){
        samtools_flags <- paste(samtools_flags, "-F2")
      } else {
        samtools_flags <- paste(samtools_flags, "-f2")
      }
    } else {
      if(keep_mapped){
        samtools_flags <- paste(samtools_flags, "-F4")
      } else {
        samtools_flags <- paste(samtools_flags, "-f4")
      }
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
      
    print(cmd)
    
    
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Z_filteraln",
      cmd = cmd,
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
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
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascet_instance=GetDefaultBascetInstance()
){
  
  #Figure out input and output file names  
  input_shards <- detect_shards_for_file(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles_R1 <- file.path(bascetRoot, input_shards) 
  
  outputFilesBAMunsorted <- make_output_shard_names(bascetRoot, outputNameBAMunsorted, "bam", num_shards)
  outputFilesBAMsorted   <- make_output_shard_names(bascetRoot, outputNameBAMsorted,   "bam", num_shards)
  
  ### What files to look for if to avoid overwriting
  outputFilesFinal <- outputFilesBAMunsorted
  if(do_sort){
    outputFilesFinal <- outputFilesBAMsorted
  }
  
  ### Verify that the input is FASTQ. check what type
  if(!is_fastq(inputFiles_R1[1])) {
    stop("Input files are not FASTQ")
  }
  
  ### Check if paired or not
  is_paired <- is_paired_fastq(inputFiles_R1[1])
  print(paste("Detect paired FASTQ:",is_paired))
  
  ### Figure out R2 names
  if(is_paired){
    inputFiles_R2 <- get_fastq_R2_from_R1(inputFiles_R1)
  } else {
    ### No R2
    inputFiles_R2 <- rep("",length(inputFiles_R1))
  }

  ### Build command: basic alignment
  final_outputFiles <- outputFilesBAMunsorted
  cmd <- c(
    shellscript_set_tempdir(bascet_instance),
    shellscript_make_bash_array("files_in_r1", inputFiles_R1),
    shellscript_make_bash_array("files_in_r2", inputFiles_R2),
    shellscript_make_bash_array("files_out_unsorted", outputFilesBAMunsorted),
    shellscript_make_bash_array("files_out_sorted", outputFilesBAMsorted),
    shellscript_make_bash_array("files_out_final", outputFilesFinal),

    ### Abort early if needed    
    if(!overwrite) helper_cancel_job_if_file_exists("${files_out_final[$TASK_ID]}"),
    
    ### For alignment
    paste(
      bascet_instance@prepend_cmd,
      "bash -c \"",
      "bwa mem", 
      useReference,
      "${files_in_r1[$TASK_ID]}",                #Align R1 FASTQ
      "${files_in_r2[$TASK_ID]}",                #Align R2 FASTQ
      "-t", numLocalThreads,
      "| ", bascet_instance@bin, "pipe-sam-add-tags",
      "| samtools view -b -o",
      "${files_out_unsorted[$TASK_ID]}",        #Each input means one output
      "\""
    )
  )
  
  
  ### Build command: sorting and indexing
  if(do_sort){
    final_outputFiles <- outputFilesBAMsorted
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
  
  
  print(cmd)
  
  if(bascet_check_overwrite_output(final_outputFiles, overwrite)) {
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Z_aln",
      cmd = cmd,
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
}









###############################################
#' Take aligned BAM file and produce Fragments.tsv.gz, compatible with Signac ATAC-seq style analysis
#' 
#' @return A job, producing a type of Fragments.tsv.gz
#' 
#' @export
BascetBam2Fragments <- function(
    bascetRoot, 
    inputName="aligned",
    outputName="fragments", 
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascet_instance=GetDefaultBascetInstance()
){
  
  input_shards <- detect_shards_for_file(bascetRoot, inputName)
  if(length(input_shards)==0){
    stop("Found no input files")
  }
  num_shards <- length(input_shards)
  inputFiles <- file.path(bascetRoot, input_shards)
  
  #One output per input pair of reads  
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "tsv.gz", num_shards)  
  
  if(bascet_check_overwrite_output(outputFiles, overwrite)) {
    RunJob(
      runner = runner, 
      jobname = "Z_bam2fragments",
      cmd = c(
        shellscript_set_tempdir(bascet_instance),
        shellscript_make_bash_array("files_in", inputFiles),
        shellscript_make_bash_array("files_out",outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) helper_cancel_job_if_file_exists("${files_out[$TASK_ID]}"),
        
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
  } else {
    new_no_job()
  }
}





###############################################
#' From aligned BAM file, compute counts per chromosome
#' 
#' @return A job, executing the counting 
#' 
#' @export
BascetCountChrom <- function(
    bascetRoot, 
    inputName="aligned",
    outputName="chromcount", 
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascet_instance=GetDefaultBascetInstance()
){
  
  input_shards <- detect_shards_for_file(bascetRoot, inputName)
  if(length(input_shards)==0){
    stop("Found no input files")
  }
  num_shards <- length(input_shards)
  inputFiles <- file.path(bascetRoot, input_shards)
  
  #One output per input alignment
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "h5", num_shards)  
  
  if(bascet_check_overwrite_output(outputFiles, overwrite)) {
    RunJob(
      runner = runner, 
      jobname = "Z_countchrom",
      cmd = c(
        shellscript_set_tempdir(bascet_instance),
        shellscript_make_bash_array("files_in", inputFiles),
        shellscript_make_bash_array("files_out",outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) helper_cancel_job_if_file_exists("${files_out[$TASK_ID]}"),
        
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
  } else {
    new_no_job()
  }
}





###############################################
#' Using Tabix, get list of sequences in a fragment file
#' 
#' TODO too raw? wrap in bascet?
#' 
#' @return The raw result of tabix
#' 
#' @export
TabixGetFragmentsSeqs <- function(
    fragpath,
    bascet_instance=GetDefaultBascetInstance()
) {
  system(paste(bascet_instance@prepend_cmd,"tabix -l ", fragpath), intern = TRUE)
}


###############################################
#' From a fragments file, get a chromatin assay for Signac  .......................... signac is dirt slow at counting! only 1 cpu used. take over the whole process
#' 
#' Case: tabix from command line uses 100% cpu only. likely not designed for large queries
#' 
#' 
#' @return A ChromatinAssay
#' 
#' @export
FragmentsToSignac <- function(
    fragpath
) {
  
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
#' 
#' @export
FragmentCountsPerChrom <- function(
    chrom_assay,
    bascet_instance=GetDefaultBascetInstance()
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
#' 
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
#' 
#' @export
ChromToSpeciesCount <- function(
    adata, 
    map_seq2strain
){
  
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
#' 
#' @export
CountGrangeFeatures <- function(
    adata, 
    grange_gene
){
  
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











################################################################################
################## SNP analysis ################################################
################################################################################

#https://academic.oup.com/bioinformatics/article/37/23/4569/6272512

#like KRAKEN, wrap cellsnp-lite

# is b needed?
# cellsnp-lite -s $BAM -b $BARCODE -O $OUT_DIR -p 10 --minMAF 0.1 --minCOUNT 100 --gzip


# singularity run ~/github/bascet/singularity/bascet.sif cellsnp-lite -s mutans_aligned.1.bam -O cellsnp.1.out -p 1 --minMAF 0.1 --minCOUNT 100 --gzip -b bclist.csv
#C1_A5_C8_A11
#F3_G6_H8_C11


# singularity run ~/mystore/bascet.sif cellsnp-lite -s mutans_aligned.1.bam -O cellsnp.1.out -p 1 --minMAF 0.1 --minCOUNT 100 --gzip -b bclist.csv
# creates output dir



###############################################
#' Align from FASTQ, generate sorted and indexed BAM file
#' 
#' TODO Should be managed by Bascet mapshard system, with automatic input conversion. unaligned file should be made temp and removed
#' 
#' @param useReference Name of the BWA reference to use
#' @param inputName xx
#' @param outputName xx
#' @param numLocalThreads Number of threads to use for each runner
#' 
#' @export
BascetRunCellSNP <- function(  
    bascetRoot, 
    numLocalThreads=1,
    inputName="aligned", 
    outputName="cellsnp", 
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascet_instance=GetDefaultBascetInstance()
){
  
  #Figure out input and output file names  
  input_shards <- detect_shards_for_file(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards) 
  
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "out", num_shards)
  
  listcellFiles <- make_output_shard_names(bascetRoot, outputName, "listcell", num_shards)
  listchromFiles <- make_output_shard_names(bascetRoot, outputName, "listchrom", num_shards)
  
  

  
  ### Build command: basic alignment
  cmd <- c(
    shellscript_set_tempdir(bascet_instance),
    shellscript_make_bash_array("files_in", inputFiles),
    shellscript_make_bash_array("files_out", outputFiles),
    shellscript_make_bash_array("files_listcell", listcellFiles),
    shellscript_make_bash_array("files_listchrom", listchromFiles),
    
    ### Abort early if needed    
    #if(!overwrite) helper_cancel_job_if_file_exists("${outputFiles[$TASK_ID]}"),  #does this work on dirs?
    if(!overwrite) helper_cancel_job_if_file_exists("${files_listcell[$TASK_ID]}"),  
    
    ### To make list of chromosome IDs
    paste(
      bascet_instance@prepend_cmd,
      "bash -c \"",
      "samtools idxstats ${files_in[$TASK_ID]} | head -n -1 | cut -f 1 > ${files_listchrom[$TASK_ID]}",
      "\""
    ),
    
    
    ### To make list of cell IDs. could avoid if we write our own cellsnp-lite
    paste(
      bascet_instance@prepend_cmd,
      "bash -c \"",
      "samtools view ${files_in[$TASK_ID]} | sed -e 's/^.*CB:Z://' | sed -e 's/\t.*//' | sort | uniq > ${files_listcell[$TASK_ID]}",
      "\""
    ),
    
    ### For SNP-calling
    paste(
      bascet_instance@prepend_cmd,

      "cellsnp-lite", 
      "-s ${files_in[$TASK_ID]}",        #Align BAM
      "-p", numLocalThreads,
      "--genotype",
      "--chrom `cat ${files_listchrom[$TASK_ID]} | paste --serial -d, - -`",
      #"--minMAF 0.1",
      #"--minCOUNT 100", #no output w these
      "--gzip",
      "-b ${files_listcell[$TASK_ID]}",
      "-O ${files_out[$TASK_ID]}"        #Each input means one output
    )
  )
  
  
  
  print(cmd)
  
  if(bascet_check_overwrite_output(outputFiles, overwrite)) {
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Z_cellsnp",
      cmd = cmd,
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
}







ReadCellSNPmatrix_one <- function(basedir) {
  
  #cellSNP.base.vcf  one line per feature
  snp_info <- read.table(file.path(basedir,"cellSNP.base.vcf.gz"))
  colnames(snp_info) <- c("chrom","pos","id","ref","alt","qual","filter","info")
  rownames(snp_info) <- sprintf(
    "%s_%s_%s_to_%s",
    snp_info$chrom, snp_info$pos, snp_info$ref, snp_info$alt
  )
  
  snp_cellid <- readLines(file.path(basedir,"cellSNP.samples.tsv"))
  length(snp_cellid)
  
  #depth of ALT alleles
  snp_ad <- Matrix::readMM(file.path(basedir,"cellSNP.tag.AD.mtx"))
  colnames(snp_ad) <- snp_cellid
  rownames(snp_ad) <- rownames(snp_info)
  
  #depth of ALT+REF alleles
  snp_dp <- Matrix::readMM(file.path(basedir,"cellSNP.tag.DP.mtx"))
  colnames(snp_dp) <- snp_cellid
  rownames(snp_dp) <- rownames(snp_info)
  
  #set feature names
  list(
    snp_ad=snp_ad,
    snp_dp=snp_dp,
    snp_info=snp_info
  )
}


###############################################
#' Read a count matrix as produced by CellSNP, but as shards
#' 
#' @return Count matrix as sparseMatrix
#' @export
ReadCellSNPmatrix <- function(
    bascetRoot, 
    inputName,
    verbose=FALSE
){
  print("Loading CellSNP file")
  
  #bascetRoot <- "/husky/henriksson/atrandi/v4_wgs_saliva1"
  #inputName <- "cellsnp"
  
  #Figure out input file names  
  input_shards <- detect_shards_for_file(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards)
  #if(tools::file_ext(inputFiles[1])!="h5"){
  #  stop("Wrong input format. should be hd5")
  #}
  
  #Load individual matrices. Sizes may not match
  list_mat <- list()
  for(f in inputFiles){
    list_one <- ReadCellSNPmatrix_one(f) ############ edited
    
    mat <- Matrix::t(list_one$snp_ad)
    
    if(verbose){
      print(dim(mat))
    }
    list_mat[[f]] <- mat
  }
  
  #Find union of features. Assume that no filtering was done such that we can assume 0
  #count for a feature not present in another shard
  all_colnames <- sort(unique(unlist(lapply(list_mat, colnames))))
  print(all_colnames)
  num_col <- length(all_colnames)
  map_name_to_i <- data.frame(row.names = all_colnames, ind=1:length(all_colnames))
  print(map_name_to_i)
  
  #Make sizes compatible
  list_resized_mat <- list()
  for(f in inputFiles){
    mat <- list_mat[[f]]
    new_mat <- MatrixExtra::emptySparse(nrow = nrow(mat), ncol = num_col, format = "R", dtype = "d")
    new_mat[1:nrow(mat), map_name_to_i[colnames(mat),]] <- MatrixExtra::as.csr.matrix(mat)  #manually look up column names!  #here, x[.,.] <- val : x being coerced from Tsparse* to CsparseMatrix
    rownames(new_mat) <- rownames(mat)
    colnames(new_mat) <- all_colnames
    # print(dim(new_mat))
    list_resized_mat[[f]] <- new_mat
  }
  
  #Concatenate matrices
  allmat <- do.call(rbind, list_resized_mat) #TODO check that above worked properly!
  
  allmat
}
