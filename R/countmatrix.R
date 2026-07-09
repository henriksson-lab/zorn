



###############################################
#' Bascet count matrix
#'
#' Sparse count matrix container used by Zorn. Rows are observations/cells and
#' columns are features.
#'
#' @param object A BascetCountMatrix object.
#' @param ... Additional arguments passed to S4 methods.
#' @slot X Sparse count matrix.
#' @slot obs Data frame with observation metadata, one row per matrix row.
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
#' Convert Bascet count storage to a Seurat-compatible counts matrix
#'
#' @noRd
BascetCountMatrixSeuratCounts <- function(mat) {
  counts <- Matrix::t(mat@X)
  counts <- methods::as(counts, "CsparseMatrix")

  if(is.null(rownames(counts))) {
    stop("Bascet count matrix has no feature names")
  }
  if(is.null(colnames(counts))) {
    stop("Bascet count matrix has no cell names")
  }
  if(anyNA(rownames(counts)) || anyDuplicated(rownames(counts))) {
    stop("Bascet count matrix feature names must be unique and non-missing")
  }
  if(anyNA(colnames(counts)) || anyDuplicated(colnames(counts))) {
    stop("Bascet count matrix cell names must be unique and non-missing")
  }

  counts
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
  SeuratObject::CreateAssayObject(counts = BascetCountMatrixSeuratCounts(mat))
}


###############################################
#' Create a Seurat object from a BascetCountMatrix
#'
#' @param counts A BascetCountMatrix object
#' @param ... Additional arguments passed to the default CreateSeuratObject method
#'
#' @return A Seurat object
#' @rawNamespace S3method(SeuratObject::CreateSeuratObject, BascetCountMatrix)
#' @export
CreateSeuratObject.BascetCountMatrix <- function(counts, ...) {
  meta.data <- counts@obs
  rownames(meta.data) <- rownames(counts@X)

  SeuratObject::CreateSeuratObject(
    counts = BascetCountMatrixSeuratCounts(counts),
    meta.data = meta.data,
    ...
  )
}


###############################################
#' Create new BascetCountMatrix
#'
#' @param X Sparse count matrix (dgCMatrix or similar)
#' @param obs Data frame of observation metadata, one row per cell
#' 
#' @return A Bascet count matrix object
#' @noRd
NewBascetCountMatrix <- function(
    X,
    obs
){
  stopifnot(nrow(obs)==nrow(X))
  if(!is.null(rownames(X))) {
    if("_index" %in% colnames(obs) &&
       !identical(as.character(obs[["_index"]]), as.character(rownames(X)))) {
      stop("Bascet count matrix obs _index is out of sync with matrix row names")
    }
    rownames(obs) <- rownames(X)
  } else {
    rownames(obs) <- NULL
  }
  new(
    "BascetCountMatrix", 
    X=X,
    obs=obs
  )
}


###############################################
# What to print to display a count matrix
#' @describeIn BascetCountMatrix Print a compact summary of a Bascet count matrix.
#' @export
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
#' @importFrom BiocGenerics rownames
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
#' @importFrom BiocGenerics `rownames<-`
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
#' @importFrom BiocGenerics colnames
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
#' @importFrom BiocGenerics `colnames<-`
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
#' @param ... Additional arguments ignored by this method
#' @param drop Ignored (kept for S4 signature compatibility)
#'
#' @return A BascetCountMatrix
#' @describeIn BascetCountMatrix Subset rows and/or columns.
#' @aliases [,BascetCountMatrix-method
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
PrepareBascetRhdf5Read <- function() {
  if(identical(Sys.getenv("HDF5_USE_FILE_LOCKING"), "")) {
    Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE")
  }
}


### Internal helper function
ReadBascetH5Shape <- function(fname, h5_index) {
  has_shape_dataset <- any(h5_index$group == "/X" & h5_index$name == "shape")
  if(has_shape_dataset) {
    # Legacy Bascet layout. Remove once pre-AnnData-style count matrices are unsupported.
    shape <- rhdf5::h5read(fname, "/X/shape")
    shape <- as.integer(shape)
    if(length(shape)==2 && all(!is.na(shape)) && all(shape >= 0)) {
      return(shape)
    }
  }

  attrs <- tryCatch(
    rhdf5::h5readAttributes(fname, "/X"),
    error = function(e) NULL
  )
  if(!is.null(attrs$shape)) {
    shape <- as.integer(attrs$shape)
    if(length(shape)==2 && all(!is.na(shape)) && all(shape >= 0)) {
      return(shape)
    }
  }

  stop("Count matrix shape could not be read from X/shape dataset or X shape attribute")
}


### Internal helper function
ReadBascetObs <- function(fname, h5_index) {
  obs_names <- h5_index$name[h5_index$group == "/obs" & h5_index$otype == "H5I_DATASET"]
  obs_values <- lapply(obs_names, function(name) rhdf5::h5read(fname, paste0("/obs/", name)))
  names(obs_values) <- obs_names
  df_obs <- as.data.frame(obs_values, optional=TRUE)
  rownames(df_obs) <- NULL

  # Older Bascet files store per-cell unclassified read counts as `_unmapped`.
  # Accept clearer aliases as well, but keep `_unmapped` for existing Zorn code.
  unclassified_aliases <- c("_unmapped", "unclassified_reads", "unmapped")
  unclassified_col <- intersect(unclassified_aliases, colnames(df_obs))[1]
  if(!is.na(unclassified_col) && !("_unmapped" %in% colnames(df_obs))) {
    df_obs$`_unmapped` <- df_obs[[unclassified_col]]
  }
  if("_unmapped" %in% colnames(df_obs) && !("unclassified_reads" %in% colnames(df_obs))) {
    df_obs$unclassified_reads <- df_obs$`_unmapped`
  }

  df_obs
}


### Internal helper function
ReadBascetCountMatrix_one <- function(
    fname,
    verbose=FALSE
){
  #check arguments
  stopifnot(file.exists(fname))
  stopifnot(is.logical(verbose))
  PrepareBascetRhdf5Read()
  
  #fname <- "/husky/henriksson/atrandi//v4_wgs_novaseq1/chromcount.1.h5"
  #fname <- "/home/mahogny/test/cnt_feature.hdf5"
  #fname <- "/husky/henriksson/atrandi//v4_wgs_novaseq1/kmer_counts.1.h5"
  
  h5_index <- rhdf5::h5ls(fname)

  indices <- rhdf5::h5read(fname, "/X/indices") + 1L
  indptr <-  rhdf5::h5read(fname, "/X/indptr")
  dat <-     rhdf5::h5read(fname, "/X/data")
  shape <-   ReadBascetH5Shape(fname, h5_index)

  mat <- Matrix::sparseMatrix(
    j=indices,
    p=indptr,
    x=dat,
    dims=shape
  )

  rownames(mat) <- as.character(rhdf5::h5read(fname, "/obs/_index"))  #names of cells
  colnames(mat) <- as.character(rhdf5::h5read(fname, "/var/_index"))  #feature names

  ### Read obs matrix
  df_obs <- ReadBascetObs(fname, h5_index)
  
  NewBascetCountMatrix(
    X=mat,
    obs=df_obs
  )
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
  bascetRoot <- normalizeBascetRoot(bascetRoot)
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
      print(dim(list_one@X))
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
#' @noRd
PrependBascetCountMatrixName <- function(
    mat,
    prependName
){
  rownames(mat@X) <- paste0(prependName, "_", rownames(mat@X))
  mat@obs$`_index` <- rownames(mat@X)
  mat
}


###############################################
#' Aggregate duplicate feature columns in a count matrix
#'
#' @noRd
AggregateDuplicateBascetFeatures <- function(mat) {
  feature_names <- as.character(colnames(mat))
  if(is.null(feature_names)) {
    stop("Bascet count matrix has no feature names")
  }
  if(anyNA(feature_names)) {
    stop("Bascet count matrix has missing feature names")
  }
  if(anyDuplicated(feature_names) == 0) {
    return(mat)
  }

  unique_features <- as.character(unique(feature_names))
  mat_summary <- Matrix::summary(mat)
  mat <- Matrix::sparseMatrix(
    i=mat_summary$i,
    j=match(feature_names[mat_summary$j], unique_features),
    x=mat_summary$x,
    dims=c(nrow(mat), length(unique_features)),
    dimnames=list(as.character(rownames(mat)), unique_features)
  )

  methods::as(mat, "CsparseMatrix")
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
  
  print("Merging count matrices")
  
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
  list_mat <- lapply(list_mat, AggregateDuplicateBascetFeatures)
  list_mat <- lapply(list_mat, function(mat) {
    rownames(mat) <- as.character(rownames(mat))
    colnames(mat) <- as.character(colnames(mat))
    mat
  })

  #Fast path for the common case: shards from one run usually have exactly the
  #same feature vector. This avoids remapping and rebuilding from triplets.
  same_colnames <- all(vapply(
    list_mat,
    function(mat) identical(colnames(mat), colnames(list_mat[[1]])),
    logical(1)
  ))

  if(same_colnames) {
    if(verbose) {
      print("Merging by direct sparse row bind")
    }
    allmat <- do.call(rbind, unname(list_mat))
    all_colnames <- as.character(sort(colnames(allmat)))
    if(!identical(colnames(allmat), all_colnames)) {
      allmat <- allmat[, all_colnames, drop=FALSE]
    }
  } else {
    #Show a progress bar
    pbar <- progress::progress_bar$new(total = numFiles)
    pbar$tick(0)

    #Find union of features
    all_colnames <- as.character(sort(unique(unlist(lapply(list_mat, colnames), use.names=FALSE))))
    if(verbose){
      print(all_colnames)
    }
    num_col <- length(all_colnames)
    col_map <- setNames(seq_along(all_colnames), all_colnames)
    if(verbose){
      print(col_map)
    }

    #General path: concatenate compressed sparse row slots directly. This avoids
    #Matrix::summary(), explicit row triplets, and a second sparse rebuild.
    list_mat <- lapply(list_mat, MatrixExtra::as.csr.matrix)
    all_p <- vector("list", length(list_mat))
    all_j <- vector("list", length(list_mat))
    all_x <- vector("list", length(list_mat))
    all_rownames <- vector("list", length(list_mat))
    nnz_offset <- 0L
    total_rows <- 0L

    for(f in seq_along(list_mat)){
      mat <- list_mat[[f]]
      col_remap <- col_map[colnames(mat)]
      remapped_j <- as.integer(unname(col_remap[mat@j + 1L]) - 1L)
      if(is.unsorted(unname(col_remap), strictly=TRUE)) {
        for(row_i in seq_len(nrow(mat))) {
          row_start <- mat@p[row_i] + 1L
          row_end <- mat@p[row_i + 1L]
          if(row_end > row_start) {
            row_ind <- row_start:row_end
            if(is.unsorted(remapped_j[row_ind], strictly=TRUE)) {
              row_order <- order(remapped_j[row_ind])
              remapped_j[row_ind] <- remapped_j[row_ind][row_order]
              mat@x[row_ind] <- mat@x[row_ind][row_order]
            }
          }
        }
      }
      all_p[[f]] <- mat@p[-1L] + nnz_offset
      all_j[[f]] <- remapped_j
      all_x[[f]] <- mat@x
      all_rownames[[f]] <- as.character(rownames(mat))
      nnz_offset <- nnz_offset + length(mat@x)
      total_rows <- total_rows + nrow(mat)
      pbar$tick()
    }

    allmat <- new(
      "dgRMatrix",
      p=as.integer(c(0L, unlist(all_p, use.names=FALSE))),
      j=as.integer(unlist(all_j, use.names=FALSE)),
      x=as.numeric(unlist(all_x, use.names=FALSE)),
      Dim=as.integer(c(total_rows, num_col)),
      Dimnames=list(as.character(unlist(all_rownames, use.names=FALSE)), all_colnames)
    )
    allmat <- methods::as(allmat, "CsparseMatrix")
  }
  
  #Concat obs
  allobs <- as.data.frame(data.table::rbindlist(list_obs, use.names=TRUE, fill=TRUE))
  rownames(allobs) <- NULL
  
  NewBascetCountMatrix(
    X=allmat,
    obs=allobs
  )
}


###############################################
#' Add metadata to a Seurat object, handling cell mismatches
#'
#' Unlike Seurat's AddMetaData, this function handles the case where
#' the metadata data.frame has fewer or more rows than cells in the object.
#' Missing cells get NA values.
#'
#' @param adata A Seurat object
#' @param metadata A data.frame with rownames matching cell names
#' @param columns Columns to add. Default: all except cell_index, taxid_index, cnt
#' @param subsetCommon If TRUE, subset the Seurat object to only cells present in both adata and metadata
#'
#' @return The Seurat object with added metadata columns
#' @export
BascetAddMetaData <- function(adata, metadata, columns=NULL, subsetCommon=FALSE) {
  if(is.null(columns)) {
    columns <- setdiff(colnames(metadata), c("cell_index","taxid_index","cnt"))
  }

  cells_in_adata <- colnames(adata)
  cells_in_meta <- rownames(metadata)

  extra_in_meta <- length(setdiff(cells_in_meta, cells_in_adata))
  extra_in_adata <- length(setdiff(cells_in_adata, cells_in_meta))

  if(subsetCommon) {
    common_cells <- intersect(cells_in_adata, cells_in_meta)
    adata <- adata[, common_cells]
    cells_in_adata <- common_cells
  } else {
    if(extra_in_meta > 0) {
      dropped <- setdiff(cells_in_meta, cells_in_adata)
      shown <- paste(utils::head(dropped, 10), collapse=", ")
      warning(paste0("BascetAddMetaData: ", extra_in_meta,
                     " cells in metadata not found in Seurat object (will be dropped): ", shown))
    }
    if(extra_in_adata > 0) {
      missing <- setdiff(cells_in_adata, cells_in_meta)
      shown <- paste(utils::head(missing, 10), collapse=", ")
      warning(paste0("BascetAddMetaData: ", extra_in_adata,
                     " cells in Seurat object not found in metadata (will be NA): ", shown))
    }
  }

  sub <- metadata[cells_in_adata, columns, drop=FALSE]

  existing_cols <- intersect(columns, colnames(adata@meta.data))
  if(length(existing_cols) > 0) {
    warning(paste0("BascetAddMetaData: replacing existing columns: ",
                   paste(existing_cols, collapse=", ")))
    adata@meta.data[, existing_cols] <- NULL
  }

  adata@meta.data <- cbind(adata@meta.data, sub)
  adata
}
