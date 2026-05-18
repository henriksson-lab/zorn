



###############################################
#' Gather all count sketches into a single count sketch matrix
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param outputName Name of output file
#' @param includeCells List of cells to process
#' @param kmerSize KMER length to use
#' @param sketchSize Size of the count sketch. Must be a power of two
#' @param overwrite Force overwriting of existing files. The default is to do nothing files exist
#' @param numThreads Number of threads to use per job. Default is the number from the runner
#' @param totalMem Total memory to allocate
#' @param runner The job manager, specifying how the command will be run (e.g. locally, or via SLURM)
#' @param bascetInstance A Bascet instance
#' 
#' @return A job to be executed, or being executed, depending on runner settings
#' @export
BascetRunCountsketch <- function( 
    bascetRoot, 
    inputName="filtered", 
    outputName="countsketch_mat.feather",  ### replace with shard!!
    includeCells=NULL,
    kmerSize=31,
    sketchSize=4096,
    overwrite=FALSE,
    numThreads=NULL,
    totalMem=NULL,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  
  #Set number of threads if not given
  if(is.null(numThreads)) {
    numThreads <- as.integer(runner@ncpu)
  }
  stopifnot(is.valid.threadcount(numThreads))
  
  #Check input arguments
  stopifnot(dir.exists(bascetRoot))
  bascetRoot <- normalizeBascetRoot(bascetRoot)
  stopifnot(is.character(outputName))
  stopifnot(is.character(inputName))
  stopifnot(is.valid.listcells(includeCells))
  stopifnot(is.numeric(kmerSize) && kmerSize == as.integer(kmerSize) && kmerSize > 0)
  stopifnot(is.numeric(sketchSize) && sketchSize == as.integer(sketchSize) && sketchSize > 0 && bitwAnd(sketchSize, sketchSize - 1) == 0)

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
    list_cell_for_shard <- rep(list(includeCells), length(inputFiles))
  }
  
  #Check memory sizes
  totalMem <- checkTotalMemArg(totalMem, runner, bascetInstance)
  
  
  if(bascetCheckOverwriteOutput(outputFile, overwrite)) {
    #Make the command
    cmd <- JobScript(
      steps = list(
        if(produce_cell_list) JobFiles("CELLFILE", list_cell_for_shard),
        JobBascetCommand(bascetInstance, list(
          "countsketch",
          if(produce_cell_list) JobArg("--cells", JobVar("CELLFILE")),
          JobMaybeArg("--memory", totalMem, format_size_bascet),
          JobArg("-@", numThreads),
          JobArg("-i", stringr::str_flatten(inputFiles, collapse = ",")),
          JobArg("--kmer-size", kmerSize),
          JobArg("--sketch-size", sketchSize),
          JobArg("-o", outputFile)
        ))
      )
    )
    
    #Run the job
    RunJob(
      runner = runner, 
      jobname = "Zcs",
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
#' @param inputName Name of countsketch matrix file. Feather and legacy CSV
#'   files are supported.
#' 
#' @return A Seurat object holding the sketch as a reduction
#' @export
BascetLoadCountSketchMatrix <- function(
    bascetRoot,
    inputName="countsketch_mat.feather"
) {
  fname <- file.path(bascetRoot, inputName)
  sketch <- ReadCountSketchMatrixFile(fname)

  cellid <- sketch$cellid
  celldepth <- sketch$celldepth
  Q <- sketch$Q
  colnames(Q) <- cellid
  rownames(Q) <- paste0("f", 1:nrow(Q))
  
  adata <- CreateSeuratObjectWithReduction(Q) #Warning: Data is of class matrix. Coercing to dgCMatrix.
  adata$celldepth <- celldepth
  
  adata
}

ReadCountSketchMatrixFile <- function(fname) {
  stopifnot(file.exists(fname))

  ext <- tolower(tools::file_ext(fname))
  if(ext %in% c("feather", "arrow", "ipc")) {
    ReadCountSketchMatrixFeather(fname)
  } else {
    ReadCountSketchMatrixCsv(fname)
  }
}

ReadCountSketchMatrixCsv <- function(fname) {
  mat <- data.table::fread(fname, header=FALSE, data.table=FALSE)
  if(ncol(mat) < 3) {
    stop("CountSketch CSV must contain cell id, depth, and at least one sketch column")
  }

  list(
    cellid = mat[[1]],
    celldepth = mat[[2]],
    Q = t(as.matrix(mat[, -(1:2), drop=FALSE]))  # Each column is one cell
  )
}

ReadCountSketchMatrixFeather <- function(fname) {
  if(!requireNamespace("arrow", quietly=TRUE)) {
    stop("Reading CountSketch Feather files requires the R package 'arrow'")
  }

  mat <- as.data.frame(arrow::read_feather(fname))
  required_cols <- c("cell_id", "depth")
  missing_cols <- setdiff(required_cols, colnames(mat))
  if(length(missing_cols) > 0) {
    stop(paste("CountSketch Feather file is missing columns:", paste(missing_cols, collapse=", ")))
  }

  sketch_cols <- grep("^cs_[0-9]+$", colnames(mat), value=TRUE)
  if(length(sketch_cols) == 0) {
    stop("CountSketch Feather file contains no cs_<index> sketch columns")
  }
  sketch_index <- as.integer(sub("^cs_", "", sketch_cols))
  sketch_cols <- sketch_cols[order(sketch_index)]

  list(
    cellid = mat[["cell_id"]],
    celldepth = mat[["depth"]],
    Q = t(as.matrix(mat[, sketch_cols, drop=FALSE]))  # Each column is one cell
  )
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
  
  m <- Matrix::Matrix(0, nrow=2, ncol=ncol(Q), sparse=TRUE)
  colnames(m) <- colnames(Q)

  pbmc <- CreateSeuratObject(counts = m, project = "a", min.cells = 0, min.features = 0)
  pbmc@reductions[[reductionName]] <- CreateDimReducObject(
    embeddings = sign(t(Q)),  ####################################### postpone??
    key = reductionName,
    assay = assay
  )
  pbmc
}


###############################################
#' Run UMAP on a count sketch reduction
#'
#' @param adata A Seurat object with a count sketch reduction
#' @param reduction Name of the reduction to use
#' @param outReduction Name of the output UMAP reduction
#' @param metric Distance metric passed to `uwot::umap`
#' @param n_neighbors Number of nearest neighbors passed to `uwot::umap`
#' @param doFast Use multi-threaded stochastic gradient updates
#' @param seed Random seed passed to `uwot::umap`
#'
#' @return Seurat object with UMAP computed
#' @export
CountSketchUMAP <- function(adata, reduction="kmersketch", outReduction="kmersketch_umap", metric="cosine", n_neighbors = 30, doFast=FALSE, seed = 42) {
  n_sgd_threads <- 0
  if(doFast) {
    print("Fast UMAP mode. Note that this means the UMAP will be different between runs (non-deterministic SGD updates)")
    n_sgd_threads <- parallel::detectCores()
  }
  
  emb <- adata@reductions[[reduction]]@cell.embeddings
  um  <- uwot::umap(
    emb, 
    metric = metric,
    n_neighbors = n_neighbors,
    nn_method = "hnsw",
    n_threads = parallel::detectCores(),
    n_sgd_threads = n_sgd_threads,
    seed = seed
  )
  colnames(um) <- c("CS_1", "CS_2")
  adata[[outReduction]] <- CreateDimReducObject(
    embeddings = um, key = "CS_", assay = DefaultAssay(adata)
  )  
  adata
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
