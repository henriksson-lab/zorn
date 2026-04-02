



###############################################
#' @export
setClass("BascetCountMatrix", slots=list(
  X="ANY", #sparse matrix
  obs="ANY"  #data frame
)
)


###############################################
#' Check that parameter is a BascetCountMatrix
#' 
#' @param x An object to test
#' @export
is.BascetCountMatrix <- function(x) {
  is(x, "BascetCountMatrix")
}


###############################################
#' Convert a BascetCountMatrix to a Seurat Assay
#'
#' @param mat A BascetCountMatrix object
#'
#' @return A Seurat AssayObject
#' @export
BascetCountMatrixToAssay <- function(mat) {
  stopifnot(is.BascetCountMatrix(mat))
  CreateAssayObject(counts = t(mat@X))
}


###############################################
#' Create new BascetCountMatrix
#'
#' @param X Sparse count matrix (dgCMatrix or similar)
#' @param obs Data frame of observation metadata, one row per cell
#' 
#' @return A Bascet count matrix object
NewBascetCountMatrix <- function(
    X,
    obs
){
  stopifnot(nrow(obs)==nrow(X))
  new(
    "BascetCountMatrix", 
    X=X,
    obs=obs
  )
}


###############################################
# What to print to display a count matrix
setMethod("show", "BascetCountMatrix", function(object) {
  cat(is(object)[[1]], "\n",
      "  Mat dim: ", nrow(object@X)," x ", ncol(object@X), "\n",
      "  Obs dim: ", nrow(object@obs), " x ", ncol(object@obs), "\n",
      sep = ""
  )
})


###############################################
#' Row names of a BascetCountMatrix (cell names)
#'
#' @param x A BascetCountMatrix object
#'
#' @return Character vector of row names
#' @export
setMethod("rownames", signature(x="BascetCountMatrix"), function(x) {
  rownames(x@X)
})


###############################################
#' Set row names of a BascetCountMatrix (cell names)
#'
#' @param x A BascetCountMatrix object
#' @param value Character vector of new row names
#'
#' @return The modified BascetCountMatrix
#' @export
setMethod("rownames<-", signature(x="BascetCountMatrix"), function(x, value) {
  rownames(x@X) <- value
  x@obs$`_index` <- value
  x
})


###############################################
#' Column names of a BascetCountMatrix (feature names)
#'
#' @param x A BascetCountMatrix object
#'
#' @return Character vector of column names
#' @export
setMethod("colnames", signature(x="BascetCountMatrix"), function(x) {
  colnames(x@X)
})


###############################################
#' Set column names of a BascetCountMatrix (feature names)
#'
#' @param x A BascetCountMatrix object
#' @param value Character vector of new column names
#'
#' @return The modified BascetCountMatrix
#' @export
setMethod("colnames<-", signature(x="BascetCountMatrix"), function(x, value) {
  colnames(x@X) <- value
  x
})


###############################################
#' Dimensions of a BascetCountMatrix
#'
#' @param x A BascetCountMatrix object
#'
#' @return Integer vector of length 2 (rows, columns) of the count matrix
#' @export
setMethod("dim", signature(x="BascetCountMatrix"), function(x) {
  dim(x@X)
})


###############################################
#' Row sums of a BascetCountMatrix (counts per cell)
#'
#' @param x A BascetCountMatrix object
#'
#' @return Named numeric vector of row sums
#' @export
setMethod("rowSums", signature(x="BascetCountMatrix"), function(x) {
  Matrix::rowSums(x@X)
})


###############################################
#' Column sums of a BascetCountMatrix (counts per feature)
#'
#' @param x A BascetCountMatrix object
#'
#' @return Named numeric vector of column sums
#' @export
setMethod("colSums", signature(x="BascetCountMatrix"), function(x) {
  Matrix::colSums(x@X)
})


###############################################
#' Subset a BascetCountMatrix
#'
#' @param x A BascetCountMatrix object
#' @param i Row (cell) indices: numeric, logical, or character vector
#' @param j Column (feature) indices: numeric, logical, or character vector
#' @param drop Ignored (kept for S4 signature compatibility)
#'
#' @return A BascetCountMatrix
#' @export
setMethod("[", signature(x="BascetCountMatrix"), function(x, i, j, ..., drop=FALSE) {
  if(!missing(i) && !missing(j)) {
    NewBascetCountMatrix(X=x@X[i, j, drop=FALSE], obs=x@obs[i, , drop=FALSE])
  } else if(!missing(i)) {
    NewBascetCountMatrix(X=x@X[i, , drop=FALSE], obs=x@obs[i, , drop=FALSE])
  } else if(!missing(j)) {
    NewBascetCountMatrix(X=x@X[, j, drop=FALSE], obs=x@obs)
  } else {
    x
  }
})



### Internal helper function
ReadBascetCountMatrix_one <- function(
    fname,
    verbose=FALSE
){
  #check arguments
  stopifnot(file.exists(fname))
  stopifnot(is.logical(verbose))
  
  #fname <- "/husky/henriksson/atrandi//v4_wgs_novaseq1/chromcount.1.h5"
  #fname <- "/home/mahogny/test/cnt_feature.hdf5"
  #fname <- "/husky/henriksson/atrandi//v4_wgs_novaseq1/kmer_counts.1.h5"
  
  h5f <- rhdf5::H5Fopen(fname)
  indices <- h5f$X$indices + 1 
  indptr <-  h5f$X$indptr
  dat <- h5f$X$data
  shape <- h5f$X$shape
  
  #shape
  #print(indices)
  
  #print(paste0("Assembling matrix, size: ", shape[1],"x",shape[2]))
  mat <- Matrix::sparseMatrix(  
    j=indices,   #i??   was j when fine
    p=indptr,
    x=dat,
    dims=h5f$X$shape
  )
  
  rownames(mat) <- h5f$obs$`_index`  #names of cells
  colnames(mat) <- h5f$var$`_index`  #feature names
  
  ### Read obs matrix
  
  #print(666)
  #print(names(h5f$obs))
  #list_obs_col <- setdiff(names(h5f$obs),"_index")
  #print(h5f$obs["_unmapped"])
  
  #df_obs <- data.frame(row.names=h5f$obs$`_index`)
  df_obs <- as.data.frame(h5f$obs)
  colnames(df_obs) <- names(h5f$obs)
  #print(head(df_obs))
  
  rhdf5::H5close()
  
  NewBascetCountMatrix(
    X=mat,
    obs=df_obs
  )
  #list(
  #  obs=df_obs,
  #  X=mat
  #)
  #mat
}


###############################################
#' Read a count matrix as produced by Bascet (hdf5 format).
#' This can be output from both BascetQueryFq and BascetCountChrom
#' 
#' @param bascetRoot The root folder where all Bascets are stored
#' @param inputName Name of input shard
#' @param verbose Print additional information, primarily to help troubleshooting
#' 
#' @return BascetCountMatrix
#' @export
ReadBascetCountMatrix <- function(
    bascetRoot, 
    inputName,
    verbose=FALSE
){
  #check arguments
  stopifnot(dir.exists(bascetRoot))
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.logical(verbose))
  
  print("Loading HDF5 files")
  
  #Figure out input file names  
  input_shards <- detectShardsForFile(bascetRoot, inputName)
  num_shards <- length(input_shards)
  if(num_shards==0){
    stop("No input files")
  }
  inputFiles <- file.path(bascetRoot, input_shards)
  if(tools::file_ext(inputFiles[1])!="h5"){
    stop("Wrong input format. should be hd5")
  }
  
  #Show a progress bar  
  pbar <- progress::progress_bar$new(total = length(inputFiles))
  pbar$tick(0)
  
  #Load individual matrices. Sizes may not match
  listInput <- list()
  #list_mat <- list()
  #list_obs <- list()
  for(f in inputFiles){
    list_one <- ReadBascetCountMatrix_one(f)
    
    if(verbose){
      print(dim(list_one$X))
    }
    listInput[[f]] <- list_one
    #list_mat[[f]] <- list_one@X
    #list_obs[[f]] <- list_one@obs
    pbar$tick()
  }
  
  MergeBascetCountMatrix(
    listInput,
    verbose=verbose
  )
}


###############################################
#' Prefix names of cells in a count matrix
#' 
PrependBascetCountMatrixName <- function(
    mat,
    prependName
){
  rownames(mat@X) <- paste0(prependName, "_", rownames(mat@X))
  mat@obs$`_index` <- rownames(mat@X)
  mat
}



###############################################
#' Merge a list of count matrices as produced by Bascet
#' 
#' @param listInput List of count matrices
#' @param prependName If TRUE, prepend name of data in list to each cell
#' @param verbose Print additional information, primarily to help troubleshooting
#' 
#' @return BascetCountMatrix
#' @export
MergeBascetCountMatrix <- function(
    listInput,
    prependName=FALSE,
    verbose=FALSE
) {
  numFiles <- length(listInput)
  
  #Prepend names to each cell before mering, if requested
  stopifnot(is.logical(prependName))
  
  if(prependName) {
    stopifnot(!is.null(names(listInput)))
    for(i in 1:numFiles) {
      prep_name <- names(listInput)[i]
      listInput[[i]] <- PrependBascetCountMatrixName(
        listInput[[i]],
        prep_name
      )
      listInput[[i]]@obs$source <- prep_name
    }
  }
  
  list_mat <- lapply(listInput, function(x) x@X)
  list_obs <- lapply(listInput, function(x) x@obs)
  
  #Show a progress bar  
  pbar <- progress::progress_bar$new(total = numFiles)
  pbar$tick(0)
  
  #Find union of features  
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
  list_resized_mat <- list()
  for(f in 1:length(list_mat)){ # inputFiles
    mat <- list_mat[[f]]
    new_mat <- MatrixExtra::emptySparse(nrow = nrow(mat), ncol = num_col, format = "C", dtype = "d")  #was format=R. but C is better, or we get warnings all the time. faster?
    new_mat[1:nrow(mat), map_name_to_i[colnames(mat),]] <- MatrixExtra::as.csr.matrix(mat)  #manually look up column n>
    rownames(new_mat) <- rownames(mat)
    colnames(new_mat) <- all_colnames
    # print(dim(new_mat))
    list_resized_mat[[f]] <- new_mat
    pbar$tick()
  }
  
  #Concatenate matrices
  allmat <- do.call(rbind, list_resized_mat) #TODO check that above worked properly!
  
  #Concat obs
  allobs <- do.call(rbind, list_obs)
  rownames(allobs) <- NULL
  
  NewBascetCountMatrix(
    X=allmat,
    obs=allobs
  )
}

