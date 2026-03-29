



###############################################
#' @export
setClass("BascetCountMatrix", slots=list(
  X="ANY", #sparse matrix
  obs="ANY"  #data frame
)
) 



###############################################
#' Create new BascetCountMatrix
#'
#' @param X The matrix
#' @param obs The observable
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




if(FALSE){
  NewBascetCountMatrix(
    mat_saliva1$X,
    mat_saliva1$obs
  )
}



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
  # 
  # 
  # #Find union of features  
  # all_colnames <- sort(unique(unlist(lapply(list_mat, colnames))))
  # if(verbose){
  #   print(all_colnames)
  # }
  # num_col <- length(all_colnames)
  # map_name_to_i <- data.frame(row.names = all_colnames, ind=1:length(all_colnames))
  # if(verbose){
  #   print(map_name_to_i)
  # }
  # 
  # #Make sizes compatible
  # list_resized_mat <- list()
  # for(f in inputFiles){
  #   mat <- list_mat[[f]]
  #   new_mat <- MatrixExtra::emptySparse(nrow = nrow(mat), ncol = num_col, format = "R", dtype = "d")
  #   new_mat[1:nrow(mat), map_name_to_i[colnames(mat),]] <- MatrixExtra::as.csr.matrix(mat)  #manually look up column names!  #here, x[.,.] <- val : x being coerced from Tsparse* to CsparseMatrix
  #   rownames(new_mat) <- rownames(mat)
  #   colnames(new_mat) <- all_colnames
  #   # print(dim(new_mat))
  #   list_resized_mat[[f]] <- new_mat
  #   pbar$tick()
  # }
  # 
  # #Concatenate matrices
  # allmat <- do.call(rbind, list_resized_mat) #TODO check that above worked properly!
  # 
  # #Concat obs
  # allobs <- do.call(rbind, list_obs)
  # 
  # list(
  #   X=allmat,
  #   obs=allobs
  # ) 
}


###############################################
#' Merge a list of count matrices as produced by Bascet
#' 
#' @param listInput List of count matrices
#' @param verbose Print additional information, primarily to help troubleshooting
#' 
#' @return BascetCountMatrix
#' @export
MergeBascetCountMatrix <- function(
    listInput,
    verbose=FALSE
){
  numFiles <- length(listInput)
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
    new_mat <- MatrixExtra::emptySparse(nrow = nrow(mat), ncol = num_col, format = "R", dtype = "d")
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
  
  NewBascetCountMatrix(
    X=allmat,
    obs=allobs
  )
  #list(
  #  X=allmat,
  #  obs=allobs
  #)
}

