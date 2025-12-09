



###############################################
#' Gather all count sketches into a single count sketch matrix
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output file
#' @param includeCells List of cells to process
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param numLocalThreads Number of threads to use per job. Default is the number from the runner
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' TODO produce a binary file format instead; gather files upon loading?
#' 
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetGatherCountSketch <- function( 
    bascetRoot, 
    inputName="countsketch", 
    outputName="countsketch_mat.csv", 
    includeCells=NULL,
    overwrite=FALSE,
    numLocalThreads=NULL,
    
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  
  #Set number of threads if not given
  if(is.null(numLocalThreads)) {
    numLocalThreads <- as.integer(runner@ncpu)
  }
  stopifnot(is.valid.threadcount(numLocalThreads))
  
  #Check input arguments
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.character(outputName))
  stopifnot(is.character(inputName))
  #includeCells todo
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  
  #Figure out input and output file names
  inputFiles <- file.path(bascetRoot, detectShardsForFile(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }
  
  outputFile <- file.path(bascetRoot, outputName)
  
  #If cell list is provided, produce a file for input (not all transform calls can handle this, so optional)
  produce_cell_list <- !is.null(includeCells)
  if(produce_cell_list) {
    #Currently using the same cell list for all shards (good idea?)
    list_cell_for_shard <- list()
    for(i in 1:length(inputFiles)){
      list_cell_for_shard[[i]] <- includeCells
    }
  }
  
  if(bascetCheckOverwriteOutput(outputFile, overwrite)) {
    #Make the command
    cmd <- c(
      #shellscript_set_tempdir(bascetInstance),
      if(produce_cell_list) shellscriptMakeFilesExpander("CELLFILE", list_cell_for_shard),
      assembleBascetCommand(bascetInstance, c(
        "countsketch",
        if(produce_cell_list) "--cells=${CELLFILE[$TASK_ID]}",
        paste0("-@=", numLocalThreads), 
        paste0("-i=", shellscriptMakeCommalist(inputFiles)),
        paste0("-o=", outputFile)
      ))
    )
    
    #Run the job
    RunJob(
      runner = runner, 
      jobname = "bascet_get_countsketch",
      bascetInstance = bascetInstance,
      cmd = cmd,
      arraysize = 1
    )  
  } else {
    new_no_job()
  }
}







###############################################
#' Compute count sketch for each cell.
#' This is a thin wrapper around BascetMapCell
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output shard
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param maxReads Max number of reads to sample per cell
#' @param kmerSize Size of the KMER to hash
#' @param sketchSize Number of dimensions to reduce to
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetComputeCountSketch___OLD <- function( 
    bascetRoot, 
    inputName="filtered", 
    outputName="countsketch", 
    overwrite=FALSE,
    maxReads=100000,  #for 5M genome, 150x2 reads, this is 6x coverage
    kmerSize=31,
    sketchSize=5000,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
) {
  BascetMapCell(
    bascetRoot=bascetRoot, 
    withfunction="_countsketch_fq", 
    inputName=inputName, 
    outputName=outputName,
    #includeCells=includeCells
    args = list(
      KMER_SIZE=format(kmerSize, scientific=FALSE),
      sketchSize=format(sketchSize, scientific=FALSE),
      MAX_READS=format(maxReads, scientific=FALSE)
    ),
    overwrite=overwrite,
    runner=runner,
    bascetInstance=bascetInstance)
}




###############################################
#' Gather all count sketches into a single count sketch matrix
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output file
#' @param includeCells List of cells to process
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param numLocalThreads Number of threads to use per job. Default is the number from the runner
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' TODO produce a binary file format instead; gather files upon loading?
#' 
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetGatherCountSketch__OLD <- function( 
    bascetRoot, 
    inputName="countsketch", 
    outputName="countsketch_mat.csv", 
    includeCells=NULL,
    overwrite=FALSE,
    
    numLocalThreads=NULL,
    
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  
  #Set number of threads if not given
  if(is.null(numLocalThreads)) {
    numLocalThreads <- as.integer(runner@ncpu)
  }
  stopifnot(is.valid.threadcount(numLocalThreads))
  
  #Check input arguments
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.character(outputName))
  stopifnot(is.character(inputName))
  #includeCells todo
  stopifnot(is.logical(overwrite))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  
  #Figure out input and output file names
  inputFiles <- file.path(bascetRoot, detectShardsForFile(bascetRoot, inputName))
  num_shards <- length(inputFiles)
  
  if(num_shards==0){
    stop("No input files")
  }
  
  outputFile <- file.path(bascetRoot, outputName)
  
  #If cell list is provided, produce a file for input (not all transform calls can handle this, so optional)
  produce_cell_list <- !is.null(includeCells)
  if(produce_cell_list) {
    #Currently using the same cell list for all shards (good idea?)
    list_cell_for_shard <- list()
    for(i in 1:length(inputFiles)){
      list_cell_for_shard[[i]] <- includeCells
    }
  }
  
  if(bascetCheckOverwriteOutput(outputFile, overwrite)) {
    #Make the command
    cmd <- c(
      #shellscript_set_tempdir(bascetInstance),
      if(produce_cell_list) shellscriptMakeFilesExpander("CELLFILE", list_cell_for_shard),
      assembleBascetCommand(bascetInstance, c(
        "countsketch",
        if(produce_cell_list) "--cells=${CELLFILE[$TASK_ID]}",
        paste0("-@=", numLocalThreads), 
        #"-@=$BASCET_TEMPDIR",
        paste0("-i=", shellscriptMakeCommalist(inputFiles)),
        paste0("-o=", outputFile)
      ))
    )
    
    #Run the job
    RunJob(
      runner = runner, 
      jobname = "bascet_get_countsketch",
      bascetInstance = bascetInstance,
      cmd = cmd,
      arraysize = 1
    )  
  } else {
    new_no_job()
  }
}








###############################################
#' Load count sketch matrix as Seurat object
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of countsketch matrix file
#' 
#' @return A Seurat object holding the sketch as a reduction
#' @export
BascetLoadCountSketchMatrix <- function(
    bascetRoot,
    inputName="countsketch_mat.csv"
) {
  fname <- file.path(bascetRoot, inputName)
  mat <- as.data.frame(data.table::fread(fname))
  
  cellid <- mat[,1]
  celldepth <- mat[,2]
  
  Q <- t(mat[,-(1:2)])  #Each column is one cell
  colnames(Q) <- cellid
  rownames(Q) <- paste0("f",1:nrow(Q))
  
  adata <- CreateSeuratObjectWithReduction(Q) #Warning: Data is of class matrix. Coercing to dgCMatrix.
  adata$celldepth <- celldepth
  
  adata
}




###############################################
#' Create a seurat object from e.g. count sketch reduction
#' 
#' @param Q Countsketch matrix
#' @param reductionName Name of reduction to store Q into
#' @param assay Name of assay to put into reduction
#' 
#' @return Seurat object holding Q as a reduction
#' @export
CreateSeuratObjectWithReduction <- function(
    Q, 
    reductionName="kmersketch", 
    assay="RNA"
){
  
  m <- matrix(nrow=2, ncol=ncol(Q))  
  colnames(m) <- colnames(Q)
  
  # as dgCMatrix 
  
  pbmc <- CreateSeuratObject(counts = m, project = "a", min.cells = 0, min.features = 0)
  pbmc@reductions[[reductionName]] <- CreateDimReducObject(
    embeddings = sign(t(Q)),  ####################################### postpone??
    key = reductionName,
    assay = assay
  )
  pbmc
}



###############################################
#' logarithmic spaced sequence; taken from emdbook library
#' 
#' @param from First value
#' @param to Last value
#' @param length.out Length of list
#' 
#' @return List of values
#' @export
lseq <- function(
    from=1, 
    to=100000, 
    length.out=6
) {
  exp(seq(log(from), log(to), length.out = length.out))
}

###############################################
#' Plot minimum number of dimensions needed to retain distance between samples
#' 
#' https://en.wikipedia.org/wiki/Johnson–Lindenstrauss_lemma
#' 
#' (1 - eps) ||u - v||^2 < ||p(u) - p(v)||^2 < (1 + eps) ||u - v||^2
#' 
#' @param listEps List of eps to plot for
#' @param minCells Plot range min cells to consider
#' @param maxCells Plot range max cells to consider
#' 
#' @return A ggplot object
#' @export
PlotJohnsonLindenstraussMinDim <- function(
    listEps, 
    minCells=10, 
    maxCells=1000000
){
  
  one_df <- data.frame(
    n_samples = lseq(minCells, maxCells, length.out=100)
  )
  
  all_df <- list()
  for(eps in listEps){
    df1 <- one_df
    df1$for_eps <- eps
    all_df[[paste(eps)]] <- df1
  }
  df <- do.call(rbind, all_df)
  
  df$n_components <- 4*log(df$n_samples) / (df$for_eps**2 / 2 - df$for_eps**3 / 3)
  
  df$eps <- paste(df$for_eps)
  
  ggplot(df, aes(n_samples, n_components, color=eps)) + 
    geom_line() +
    scale_x_log10() +
    scale_y_log10() +
    xlab("Number of cells")+
    ylab("Dimensions") +
    theme_bw()
}














if(FALSE){
  library(future)
  plan("multicore", workers = 10)
  
  pbmc <- CreateSeuratObjectWithReduction(Q[,1:1000])
  pbmc
  reductionName <- "kmersketch"
  pbmc <- RunUMAP(
    pbmc, 
    dims = 1:ncol(pbmc@reductions[[reductionName]]@cell.embeddings), 
    reduction = reductionName,
    metric = "cosine"
  )  
  #todo cosine distance??
  ##pbmc <- FindNeighbors(pbmc, reduction = "kmersketch", annoy.metric = "cosine")
}