 

################################################################################
################## Alignment ###################################################
################################################################################


###############################################
#' Figure out if a BAM-file is a paired alignment or not
#' 
#' @param fname Name of BAM-file
#' 
#' @export
isBamPairedAlignment <- function(
    fname,
    bascetInstance=GetDefaultBascetInstance()
){

  #example
  #@PG     ID:bwa  PN:bwa  VN:0.7.18-r1243-dirty   CL:bwa mem /home/m/mahogny/mystore/atrandi/bwa_ref/combo_hx/all.fa.gz /pfs/proj/nobackup/fs/projnb10/hpc2nstor2024-027/atrandi/v2_wgs_saliva1/asfq.1.R1.fq.gz /pfs/proj/nobackup/fs/projnb10/hpc2nstor2024-027/atrandi/v2_wgs_saliva1/asfq.1.R2.fq.gz -t 10
  
  #Note: can also see this from the BAM header, but a real mess to parse out!
  # samtools view unsorted_aligned.1.bam  -H
  # @PG	ID:bwa	PN:bwa	VN:0.7.17-r1198-dirty	CL:bwa mem /husky/fromsequencer/241210_joram_rnaseq/ref/all.fa /husky/henriksson/atrandi/v2_rnaseq5//asfq.1.R1.fq.gz /husky/henriksson/atrandi/v2_rnaseq5/asfq.1.R2.fq.gz -t 20
  
  cmd <- paste(
    bascetInstance@prependCmd, 
    "samtools view",
    fname,
    "-h | head -n 10000 | ", #100 was not enough to be reliable
    bascetInstance@prependCmd," samtools view -c -f 1"
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
    bascetInstance=GetDefaultBascetInstance()
){
  
  if(!file.exists(genomeFile)){
    stop("Could not find genome FASTA file")
  }
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Z_bwa_index",
      bascetInstance = bascetInstance,
      cmd = c(
        ### For sorting
        paste(
          bascetInstance@prependCmd,
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
    bascetInstance=GetDefaultBascetInstance()
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
    bascetInstance = bascetInstance,
    cmd = c(
      paste(
        "mdkir -p", outDir
      ),
      paste(
        bascetInstance@prependCmd,
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
isFastq <- function(fname) {
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
isPairedFastq <- function(fname) {
  if(isFastq(fname)){
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
getFastqR2fromR1 <- function(fname) {

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
    bascetInstance=GetDefaultBascetInstance()
){
  
  #Figure out input and output file names
  input_shards <- detectShardsForFile(bascetRoot, inputName) 
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards) 
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "bigwig", num_shards)
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    
    ### Build command
    cmd <- c(
      #shellscript_set_tempdir(bascetInstance),
      shellscriptMakeBashArray("files_in", inputFiles),
      shellscriptMakeBashArray("files_out", outputFiles),
      
      ### For sorting
      paste(
        bascetInstance@prependCmd,
        "bamCoverage",
        "-b ${files_in[$TASK_ID]}",     #Each job takes a single output
        "-o ${files_out[$TASK_ID]}"     #Each job produces a single output
      )
    )
    
    
    
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Z_tobw",
      bascetInstance = bascetInstance,
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
#' @param keepMapped Keep the mapped reads (TRUE) or unmapped (FALSE)
#' 
#' @export
BascetFilterAlignment <- function(  
    bascetRoot, 
    numLocalThreads=1,
    inputName, 
    outputName,
    keepMapped=FALSE,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  
  
  #Figure out input and output file names
  input_shards <- detectShardsForFile(bascetRoot, inputName) 
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards) 
  
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "bam", num_shards)

  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    
    ### Set flags for samtools. Depends on if paired alignment or not
    samtools_flags <- ""
    is_paired_al <- isBamPairedAlignment(inputFiles[1], bascetInstance)
    print(paste("Detect paired alignment: ",is_paired_al))
    
    if(is_paired_al) {
      if(keepMapped){
        samtools_flags <- paste(samtools_flags, "-F2")
      } else {
        samtools_flags <- paste(samtools_flags, "-f2")
      }
    } else {
      if(keepMapped){
        samtools_flags <- paste(samtools_flags, "-F4")
      } else {
        samtools_flags <- paste(samtools_flags, "-f4")
      }
    }
    
    ### Build command
    cmd <- c(
      #shellscript_set_tempdir(bascetInstance),
      shellscriptMakeBashArray("files_in", inputFiles),
      shellscriptMakeBashArray("files_out", outputFiles),
      
      ### For sorting
      paste(
        bascetInstance@prependCmd,
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
      bascetInstance = bascetInstance,
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
    aligner="BWA",
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  
  #Figure out input and output file names  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles_R1 <- file.path(bascetRoot, input_shards) 
  
  outputFilesBAMunsorted  <- makeOutputShardNames(bascetRoot, outputNameBAMunsorted, "bam", num_shards)
  outputFilesBAMsorted    <- makeOutputShardNames(bascetRoot, outputNameBAMsorted,   "bam", num_shards)
  outputFilesBAMsortedPRE <- stringr::str_remove(outputFilesBAMsorted, "\\.bam$")
  
  ### What files to look for if to avoid overwriting
  outputFilesFinal <- outputFilesBAMunsorted
  if(do_sort){
    outputFilesFinal <- outputFilesBAMsorted
  }
  
  ### Verify that the input is FASTQ. check what type
  if(!isFastq(inputFiles_R1[1])) {
    stop("Input files are not FASTQ")
  }
  
  ### Check if paired or not
  is_paired <- isPairedFastq(inputFiles_R1[1])
  print(paste("Detect paired FASTQ:",is_paired))
  
  ### Figure out R2 names
  if(is_paired){
    inputFiles_R2 <- getFastqR2fromR1(inputFiles_R1)
  } else {
    ### No R2
    inputFiles_R2 <- rep("",length(inputFiles_R1))
  }

  
  if(aligner=="BWA"){
    
    cmd_align <- paste(
      bascetInstance@prependCmd,
      "bash -c \"",
      "bwa mem", 
      useReference,
      "${files_in_r1[$TASK_ID]}",                #Align R1 FASTQ
      "${files_in_r2[$TASK_ID]}",                #Align R2 FASTQ
      "-t", numLocalThreads,
      "| ", bascetInstance@bin, "pipe-sam-add-tags",
      "| samtools view - -S -b -o",
      "${files_out_unsorted[$TASK_ID]}",        #Each input means one output
      "\""
    )
    
    print(cmd_align)
    
  } else if(aligner=="STAR") {
    

    cmd_mkdir_starlog <- paste(  ########## TODO delete all of this directory?
      paste("mkdir -p STARlog")
    )
    
    #### STAR needs its own tempdir. Make sure they don't collide.
    #TODO Should we rather put in subdir of user request?
    star_temp_dir <- "STARlog/_STARtmp.${TASK_ID}"
    
    #ensure any old STAR temp directory is removed.
    #the naming should make this command safe
    cmd_delete_temp <- paste(
      paste("rm -Rf",star_temp_dir)
    )
    
    cmd_star <- paste(
      bascetInstance@prependCmd,
      "bash -c \"",
      "STAR", 
      "--genomeDir", useReference,
      "--readFilesIn ${files_in_r1[$TASK_ID]} ${files_in_r2[$TASK_ID]}",  #Align R1 and R2 FASTQ
      "--runThreadN", numLocalThreads,
      "--outSAMtype SAM",  #this implies Unsorted
      "--outSAMunmapped Within",
      "--outSAMattributes Standard",
      paste("--outTmpDir",star_temp_dir), 
      "--outFileNamePrefix ./STARlog/${TASK_ID}_",
      "--outStd SAM",
      "--readFilesCommand zcat",
      "| ", bascetInstance@bin, "pipe-sam-add-tags",
      "| samtools view - -S -b -o",
      "${files_out_unsorted[$TASK_ID]}",        #Each input means one output
      "\""
    )
    
    cmd_align <- c(
      cmd_mkdir_starlog,
      cmd_delete_temp,
      cmd_star
    )
  } else {
    stop(paste("Unknown aligner: ", aligner))
  }
  
  
  ### Build command: basic alignment
  final_outputFiles <- outputFilesBAMunsorted
  cmd <- c(
    #shellscript_set_tempdir(bascetInstance),
    shellscriptMakeBashArray("files_in_r1", inputFiles_R1),
    shellscriptMakeBashArray("files_in_r2", inputFiles_R2),
    shellscriptMakeBashArray("files_out_unsorted", outputFilesBAMunsorted),
    shellscriptMakeBashArray("files_out_sorted", outputFilesBAMsorted),
    shellscriptMakeBashArray("files_out_sorted_pre", outputFilesBAMsortedPRE),
    shellscriptMakeBashArray("files_out_final", outputFilesFinal),

    ### Abort early if needed    
    if(!overwrite) shellscriptCancelJobIfFileExists("${files_out_final[$TASK_ID]}"),
    
    ### For alignment
    "echo Aligning ${files_in_r1[$TASK_ID]} + ${files_in_r2[$TASK_ID]} ==> ${files_out_unsorted[$TASK_ID]}",
    cmd_align
  )
  
  
  ### Build command: sorting and indexing
  if(do_sort){
    final_outputFiles <- outputFilesBAMsorted
    cmd <- c(
      cmd,
      
      ### For sorting
      "echo Sorting ${files_out_unsorted[$TASK_ID]}  ==>  ${files_out_sorted[$TASK_ID]}",
      paste(
        bascetInstance@prependCmd,
        "samtools sort", 
#        "-f",                                 #Output file name is given as a whole (not prefix)
        "-@",numLocalThreads,                 #Number of threads to use
        "${files_out_unsorted[$TASK_ID]}",    #Each job produces a single output
        "${files_out_sorted_pre[$TASK_ID]}"       #Each job produces a single output
#        "${files_out_sorted[$TASK_ID]}"       #Each job produces a single output
      ),
      
      ### For indexing
      "echo Indexing ${files_out_sorted[$TASK_ID]}",
      paste(
        bascetInstance@prependCmd,
        "samtools index", 
        "${files_out_sorted[$TASK_ID]}"        #Each job produces a single output
      )      
    )
  }
  
  
  print(cmd)
  
  if(bascetCheckOverwriteOutput(final_outputFiles, overwrite)) {
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Z_aln",
      bascetInstance = bascetInstance,
      cmd = cmd,
      arraysize = num_shards
    )
  } else {
    new_no_job()
  }
}


#STAR on parse Could not parse UMI from read name

#
# singularity run ~/mystore//bascet.sif   bash -c " STAR --genomeDir /home/m/mahogny/mystore/dataset/ref_human_virusas/GRCh38 --readFilesIn ${files_in_r1[$TASK_ID]} ${files_in_r2[$TASK_ID]} --runThreadN 10 --outSAMtype SAM --outSAMunmapped Within --outSAMattributes Standard |  bascet pipe-sam-add-tags | samtools view -b -o ${files_out_unsorted[$TASK_ID]} "
# singularity run ~/mystore//bascet.sif   bash -c " STAR --genomeDir /home/m/mahogny/mystore/dataset/ref_human_virusas/GRCh38 --readFilesIn fastp.1.R1.fq.gz fastp.1.R2.fq.gz --runThreadN 10 --outSAMtype SAM --outSAMunmapped Within --outSAMattributes Standard --outStd SAM > foo.sam "















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
    bascetInstance=GetDefaultBascetInstance()
){
  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  if(length(input_shards)==0){
    stop("Found no input files")
  }
  num_shards <- length(input_shards)
  inputFiles <- file.path(bascetRoot, input_shards)
  
  #One output per input pair of reads  
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "tsv.gz", num_shards)  
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    RunJob(
      runner = runner, 
      jobname = "Z_bam2fragments",
      bascetInstance = bascetInstance,
      cmd = c(
        #shellscript_set_tempdir(bascetInstance),
        shellscriptMakeBashArray("files_in", inputFiles),
        shellscriptMakeBashArray("files_out",outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),
        
        paste(
          bascetInstance@prependCmd,
          bascetInstance@bin, 
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
    minMatching=0,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  if(length(input_shards)==0){
    stop("Found no input files")
  }
  num_shards <- length(input_shards)
  inputFiles <- file.path(bascetRoot, input_shards)
  
  #One output per input alignment
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "h5", num_shards)  
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    RunJob(
      runner = runner, 
      jobname = "Z_countchrom",
      bascetInstance = bascetInstance,
      cmd = c(
        #shellscript_set_tempdir(bascetInstance),
        shellscriptMakeBashArray("files_in", inputFiles),
        shellscriptMakeBashArray("files_out",outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),
        
        paste(
          bascetInstance@prependCmd,
          bascetInstance@bin, 
          "countchrom",
          "--min-matching ",minMatching,
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
    bascetInstance=GetDefaultBascetInstance()
) {
  system(paste(bascetInstance@prependCmd,"tabix -l ", fragpath), intern = TRUE)
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
    bascetInstance=GetDefaultBascetInstance()
){
  
  #Figure out where the fragment file is
  fr <- Fragments(chrom_assay)[[1]]
  #fr@cells ## also possible!
  
  #Get the name of chromosomes
  all_seqid <- TabixGetFragmentsSeqs(fr@path, bascetInstance)
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
  #includeCells <- sort(intersect(colnames(adata), colnames(chrom_assay)))
  #chrom_assay <- chrom_assay[,includeCells]
  #adata <- adata[,includeCells]
  
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
    mapSeq2strain
){
  
  mat_cnt <- adata@assays[[DefaultAssay(adata)]]$counts
  
  unique_strains <- unique(mapSeq2strain$strain)
  strain_cnt <- matrix(NA, ncol=ncol(mat_cnt), nrow=length(unique_strains))
  rownames(strain_cnt) <- unique_strains
  colnames(strain_cnt) <- colnames(mat_cnt)
  for(i in 1:nrow(strain_cnt)){
    cur_cols <- mapSeq2strain$id[mapSeq2strain$strain %in% rownames(strain_cnt)[i]]
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
#' @param grangeGene A grange object telling where to count
#' @return A Seurat AssayObject
#' 
#' @export
CountGrangeFeatures <- function(
    adata, 
    grangeGene
){
  
  gene_counts <- FeatureMatrix(
    fragments = Fragments(adata),
    features = grangeGene,
    cells = colnames(adata)
  )
  rownames(gene_counts) <- grangeGene$Name ## paste(grangeGene$gene_biotype,grangeGene$Name) #Name in an easy manner
  
  # create a new assay using the MACS2 peak set 
  CreateAssayObject(
    counts = gene_counts[!duplicated(rownames(gene_counts)),],
    min.cells = 0,  #filter later!!!
    min.features = 0
  )
}













###############################################
#' From aligned BAM file, compute counts per feature
#' 
#' 
#' 
#' @return A job, executing the counting 
#' 
#' @export
BascetCountFeature <- function(
    bascetRoot,
    inputName="aligned",
    outputName="featurecount", 
    gffFile,
    useFeature="gene",
    attrGeneId="gene_id",  #"ID" for yersinia
    attrGeneName="name",
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(), 
    bascetInstance=GetDefaultBascetInstance()
){
  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  if(length(input_shards)==0){
    stop("Found no input files")
  }
  num_shards <- length(input_shards)
  inputFiles <- file.path(bascetRoot, input_shards)
  
  #One output per input alignment
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "h5", num_shards)  
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    RunJob(
      runner = runner, 
      jobname = "Z_countf",
      bascetInstance = bascetInstance,
      cmd = c(
        shellscriptMakeBashArray("files_in", inputFiles),
        shellscriptMakeBashArray("files_out",outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),
        
        paste(
          bascetInstance@prependCmd,
          bascetInstance@bin, 
          "countfeature",
          paste("-g", gffFile),
          paste("--use-feature", useFeature),
          paste("--attr-id", attrGeneId),
          paste("--attr-name", attrGeneName),
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
    bascetInstance=GetDefaultBascetInstance()
){
  
  #Figure out input and output file names  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards) 
  
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "out", num_shards)
  
  listcellFiles <- makeOutputShardNames(bascetRoot, outputName, "listcell", num_shards)
  listchromFiles <- makeOutputShardNames(bascetRoot, outputName, "listchrom", num_shards)
  
  

  
  ### Build command: basic alignment
  cmd <- c(
    #shellscript_set_tempdir(bascetInstance),
    shellscriptMakeBashArray("files_in", inputFiles),
    shellscriptMakeBashArray("files_out", outputFiles),
    shellscriptMakeBashArray("files_listcell", listcellFiles),
    shellscriptMakeBashArray("files_listchrom", listchromFiles),
    
    ### Abort early if needed    
    #if(!overwrite) shellscriptCancelJobIfFileExists("${outputFiles[$TASK_ID]}"),  #does this work on dirs?
    if(!overwrite) shellscriptCancelJobIfFileExists("${files_listcell[$TASK_ID]}"),  
    
    ### To make list of chromosome IDs
    paste(
      bascetInstance@prependCmd,
      "bash -c \"",
      "samtools idxstats ${files_in[$TASK_ID]} | head -n -1 | cut -f 1 > ${files_listchrom[$TASK_ID]}",
      "\""
    ),
    
    
    ### To make list of cell IDs. could avoid if we write our own cellsnp-lite
    paste(
      bascetInstance@prependCmd,
      "bash -c \"",
      "samtools view ${files_in[$TASK_ID]} | sed -e 's/^.*CB:Z://' | sed -e 's/\t.*//' | sort | uniq > ${files_listcell[$TASK_ID]}",
      "\""
    ),
    
    ### For SNP-calling
    paste(
      bascetInstance@prependCmd,

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
  
  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    #Produce the script and run the job
    RunJob(
      runner = runner, 
      jobname = "Z_cellsnp",
      bascetInstance = bascetInstance,
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
    listCells=NULL,
    verbose=FALSE
){
  print("Loading CellSNP file")
  
  #bascetRoot <- "/husky/henriksson/atrandi/v4_wgs_saliva1"
  #inputName <- "cellsnp"
  
  #Figure out input file names  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards)
  #if(tools::file_ext(inputFiles[1])!="h5"){
  #  stop("Wrong input format. should be hd5")
  #}
  
  #Load individual matrices. Sizes may not match
  if(verbose){
    print("Loading matrices")
  }
  list_mat <- list()
  for(f in inputFiles){
    list_one <- ReadCellSNPmatrix_one(f) 
    
    mat <- Matrix::t(list_one$snp_ad)
    
    if(!is.null(listCells)){
      rownames(mat)
      mat <- mat[rownames(mat) %in% listCells,,drop=FALSE]
    }
    
    if(verbose){
      print(dim(mat))
    }
    list_mat[[f]] <- mat
  }
  
  #Find union of features. Assume that no filtering was done such that we can assume 0
  #count for a feature not present in another shard
  if(verbose){
    print("Union of features across matrices")
  }
  all_colnames <- sort(unique(unlist(lapply(list_mat, colnames))))
  if(verbose){
    print(all_colnames)
  }
  num_col <- length(all_colnames)
  map_name_to_i <- data.frame(row.names = all_colnames, ind=1:length(all_colnames))
  if(verbose){
    print(map_name_to_i)
  }
  
  #Make sizes compatible
  if(verbose){
    print("Concatenating")
  }
  list_resized_mat <- list()
  for(f in inputFiles){
    if(verbose){
      print(f)
    }
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
