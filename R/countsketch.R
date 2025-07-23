

###############################################
#' Compute count sketch for each cell.
#' This is a thin wrapper around BascetMapCell
#' 
#' @inheritParams template_BascetFunction
#' @return TODO
#' @export
BascetComputeCountSketch <- function( 
    bascetRoot, 
    inputName="filtered", 
    outputName="countsketch", 
    #includeCells=NULL,
    overwrite=FALSE,
    maxReads=100000,  #for 5M genome, 150x2 reads, this is 6x coverage
    kmerSize=31,
    sketch_size=5000,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  BascetMapCell(
    bascetRoot=bascetRoot, 
    withfunction="_countsketch_fq", 
    inputName=inputName, 
    outputName=outputName,
    #includeCells=includeCells
    args = list(
      KMER_SIZE=format(kmerSize, scientific=FALSE),
      SKETCH_SIZE=format(sketch_size, scientific=FALSE),
      MAX_READS=format(maxReads, scientific=FALSE)
    ),
    overwrite=overwrite,
    runner=runner,
    bascetInstance=bascetInstance)
}




###############################################
#' Gather all count sketches into a single count sketch matrix
#' 
#' TODO binary file format
#' 
#' @inheritParams template_BascetFunction
#' @return A job
#' @export
BascetGatherCountSketch <- function( 
    bascetRoot, 
    inputName="countsketch", 
    outputName="countsketch_mat.csv", 
    includeCells=NULL,
    overwrite=FALSE,
    runner=GetDefaultBascetRunner(),
    bascetInstance=GetDefaultBascetInstance()
){
  
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
      paste(
        bascetInstance@prependCmd,
        bascetInstance@bin, 
        "countsketch",
        if(produce_cell_list) "--cells ${CELLFILE[$TASK_ID]}",
        "-t $BASCET_TEMPDIR",
        "-i", shellscriptMakeCommalist(inputFiles),
        "-o", outputFile
      )
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
#' @return A seurat object
#' @export
BascetLoadCountSketchMatrix <- function(
    bascetRoot,
    inputName="countsketch_mat.csv"
) {
  #"/husky/henriksson/atrandi/v2_wgs_novaseq1/countsketch_mat.csv"
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
#' @return Seurat object; make an assay?
#' @export
CreateSeuratObjectWithReduction <- function(Q, reduction_name="kmersketch", assay="RNA"){
  
  m <- matrix(nrow=2, ncol=ncol(Q))  
  colnames(m) <- colnames(Q)
  
  # as dgCMatrix 
  
  pbmc <- CreateSeuratObject(counts = m, project = "a", min.cells = 0, min.features = 0)
  pbmc@reductions[[reduction_name]] <- CreateDimReducObject(
    embeddings = sign(t(Q)),  ####################################### postpone??
    key = reduction_name,
    assay = assay
  )
  pbmc
}



###############################################
#' logarithmic spaced sequence, from emdbook library
lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

###############################################
#' Plot minimum number of dimensions needed to retain distance between samples
#' 
#' https://en.wikipedia.org/wiki/Johnson–Lindenstrauss_lemma
#' 
#' (1 - eps) ||u - v||^2 < ||p(u) - p(v)||^2 < (1 + eps) ||u - v||^2
#' 
#' @param list_eps List of eps to plot for
#' @export
PlotJohnsonLindenstraussMinDim <- function(list_eps, min_cells=10, max_cells=1000000){
  
  one_df <- data.frame(
    n_samples = lseq(min_cells, max_cells, length.out=100)
  )
  
  all_df <- list()
  for(eps in list_eps){
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
  reduction_name <- "kmersketch"
  pbmc <- RunUMAP(
    pbmc, 
    dims = 1:ncol(pbmc@reductions[[reduction_name]]@cell.embeddings), 
    reduction = reduction_name,
    metric = "cosine"
  )  
  
  
  #todo cosine distance??
  
  ##pbmc <- FindNeighbors(pbmc, reduction = "kmersketch", annoy.metric = "cosine")
  
}
