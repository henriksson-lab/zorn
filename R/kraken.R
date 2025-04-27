

###############################################
# TODO finish properly
SpeciesCorrMatrix <- function(
    adata
){
  cnt <- adata@assays[[DefaultAssay(adata)]]$counts
  cnt <- cnt[rowSums(cnt)>0,]
  
  list_species <- rownames(cnt)
  print(list_species)
  
  all_comp <- NULL
  for(i in seq_along(list_species)){
    for(j in seq_along(list_species)){
      #print(paste(i,j))
      if(i==j) {
        
      } else {
        
        df <- data.frame(
          x=factor(cnt[i,]>0, levels=c("TRUE","FALSE")),
          y=factor(cnt[j,]>0, levels=c("TRUE","FALSE"))
        )
        
        df <- data.frame(
          x=factor(cnt[i,]>median(cnt[i,]), levels=c("TRUE","FALSE")),
          y=factor(cnt[j,]>median(cnt[j,]), levels=c("TRUE","FALSE"))
        )
        
        ft <- fisher.test(table(df))
        
        all_comp <- rbind(all_comp,
                          data.frame(
                            i=list_species[i], 
                            j=list_species[j], 
                            p=ft$p.value
                          ))
      }
    }
  }
  ggplot(all_comp, aes(i,j,fill = -log(p))) + geom_tile() + theme_bw()
  #all_comp
  #egg::ggarrange(plots = all_plots, nrow = nrow(cnt))  
}






###############################################
#' Run KRAKEN2 for each cell
#'
#' @param useKrakenDB Path to KRAKEN2 database
#' @export
#' 
BascetRunKraken <- function(
    bascetRoot,
    useKrakenDB="/data/henlab/kraken/standard-8",
    numLocalThreads=1,
    inputName="asfq", ######### should be able to take filtered and pipe to kraken if needed  "filtered"
    outputName="kraken_out",
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
  
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "kraken_out", num_shards) 
  
  ### Check if paired or not
  is_paired <- is_paired_fastq(inputFiles_R1[1])
  print(paste("Detect paired FASTQ:",is_paired))
  
  ### Figure out R2 names
  if(is_paired){
    inputFiles_R2 <- get_fastq_R2_from_R1(inputFiles_R1)
  }
  
  if(bascet_check_overwrite_output(outputFiles, overwrite)) {
    #Run the job
    RunJob(
      runner = runner, 
      jobname = paste0("Z_kraken_fq"),
      cmd = c(
        shellscript_set_tempdir(bascet_instance),
        shellscript_make_bash_array("files_in_R1",inputFiles_R1),
        if(is_paired) shellscript_make_bash_array("files_in_R2",inputFiles_R2),
        shellscript_make_bash_array("files_out",outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) helper_cancel_job_if_file_exists("${files_out[$TASK_ID]}"),
        
        paste(
          bascet_instance@prepend_cmd,
          "kraken2",
          "--db", useKrakenDB,
          "--threads", numLocalThreads,
          "--output ${files_out[$TASK_ID]}",
          if(is_paired) "--paired",
          "${files_in_R1[$TASK_ID]}",
          if(is_paired) "${files_in_R2[$TASK_ID]}"
        )
      ),
      arraysize = num_shards
    )  
  } else {
    new_no_job()
  }
}




###############################################
#' Produce a count matrix of taxonomy IDs from KRAKEN output
#' 
#' TODO: should run kraken in the same run
#' 
#' @inheritParams template_BascetFunction
#' @param useKrakenDB description
#' @param numLocalThreads description
#' @param inputName description
#' @param outputName description
#' @export
BascetMakeKrakenCountMatrix <- function(
    bascetRoot,
    numLocalThreads=1,
    inputName="kraken_out", ######### should be able to take filtered and pipe to bwa if needed  "filtered"
    outputName="kraken", 
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
  
  outputFiles <- make_output_shard_names(bascetRoot, outputName, "h5", num_shards) 
  
  if(bascet_check_overwrite_output(outputFiles, overwrite)) {
    #Run the job
    RunJob(
      runner = runner, 
      jobname = paste0("Z_kraken_mat"),
      cmd = c(
        shellscript_set_tempdir(bascet_instance),
        shellscript_make_bash_array("files_in",inputFiles),
        shellscript_make_bash_array("files_out",outputFiles),
        
        ### Abort early if needed    
        if(!overwrite) helper_cancel_job_if_file_exists("${files_out[$TASK_ID]}"),

        paste(
          bascet_instance@prepend_cmd,
          bascet_instance@bin, 
          "kraken",
          "-t $BASCET_TEMPDIR",
          "-i ${files_in[$TASK_ID]}",
          "-o ${files_out[$TASK_ID]}"
          )
      ),
      arraysize = num_shards
    )  
  } else {
    new_no_job()
  }
}




#' ###############################################
#' #' Read a KRAKEN2 count matrix as produced by Bascet (hdf5 format)
#' #' 
#' #' @param fname Full name of the HDF5 count matrix file
#' #' @return Counts as a sparseMatrix
#' #' @export
#' ReadBascetKrakenMatrix <- function(
#'     bascetRoot,
#'     inputName="kraken"
#' ){
#'   #Figure out input file names  
#'   input_shards <- detect_shards_for_file(bascetRoot, inputName)
#'   num_shards <- length(input_shards)
#'   if(num_shards==0){
#'     stop("No input files")
#'   }
#'   inputFiles <- file.path(bascetRoot, input_shards)
#'   
#'   #Load individual matrices. these may not have compatible sizes
#'   list_mat <- list()
#'   for(f in inputFiles){
#'     mat <- ReadBascetKrakenMatrix_one(f)
#'     list_mat[[f]] <- mat
#'   }
#'   
#'   ## Resize the matrices to have matching number of rows.
#'   ## Could also assemble them at this step if faster (probably not)
#'   max_row <- max(sapply(list_mat, nrow))
#'   list_resized_mat <- list()
#'   for(f in inputFiles){
#'     mat <- list_mat[[f]]
#'     new_mat <- MatrixExtra::emptySparse(nrow = max_row, ncol = ncol(mat), format = "R", dtype = "d")
#'     new_mat[1:nrow(mat), 1:ncol(mat)] <- mat
#'     rownames(new_mat) <- rownames(mat)
#'     colnames(new_mat) <- colnames(mat)
#'     list_resized_mat[[f]] <- new_mat
#'   }
#' 
#'   all_mat <- do.call(cbind, list_resized_mat)
#'   all_mat
#' }
#' 
#' 
#' ###############################################
#' #' Internal function, to load a single kraken matrix
#' ReadBascetKrakenMatrix_one <- function(
#'     fname
#' ){
#'   h5f <- rhdf5::H5Fopen(fname)
#'   indices <- h5f$X$indices+1
#'   indptr <-  h5f$X$indptr
#'   dat <- h5f$X$data
#'   shape <- h5f$X$shape
#'   
#'   #print(paste0("Assembling matrix, size: ", shape[1],"x",shape[2]))
#'   mat <- Matrix::sparseMatrix(  
#'     j=indices, 
#'     p=indptr,
#'     x=dat,
#'     dims=h5f$X$shape
#'   )
#'   
#'   rownames(mat) <- h5f$obs$`_index`
#'   
#'   rhdf5::H5close()
#'   
#'   mat <- Matrix::t(mat)
#'   ##Note that taxid 0 is added, but with index 1 in R. Thus need to remove the first column
#'   unident <- mat[-1,]
#'   mat <- mat[-1,]
#'   #  dim(mat)
#' }





###############################################
#' For a KRAKEN2 count matrix, return consensus taxID for each cell as metadata
#' 
#' @param mat A count matrix, typically in sparse format
#' @return A data.frame holding cellID and consensus taxID
#' @export
KrakenFindConsensusTaxonomy <- function(
    mat
){
  
  #taxid 135 9167  7942
  
  
  #turn into triplet representation
  library(Matrix)
  #M <- Matrix::Matrix(mat, sparse = TRUE)
  M <- as(mat, "TsparseMatrix")
  
  #Note that indices are zero-based within Matrix package
  df <- data.frame(
    cell_index = M@i + 1,
    taxid_index = M@j + 1,
    cnt = M@x
  )
  
  col_taxid <- stringr::str_remove(colnames(M),"taxid_")
  df$taxid <- col_taxid[df$taxid_index]
  
  #Not all taxid's will map. so we need to look them up first, then discard some of them
  taxonomizr::prepareDatabase(getAccessions=FALSE)
  
  for_taxid <- unique(df$taxid)
  df_taxid <- data.frame(taxonomizr::getTaxonomy(
    for_taxid,
    desiredTaxa = c(
      "phylum", "class", "order", "family", "genus","species")
  ))
  df_taxid$taxid <- for_taxid
  df_taxid <- df_taxid[!is.na(df_taxid$species),] #keep entries with known species
  
  #Get taxid with most counts per cell. Only keep those that are for species. Note that more than one taxid can be reported!
  df <- df[df$taxid %in% df_taxid$taxid,]
  max_taxid_per_cell <- merge(df,sqldf::sqldf("select max(cnt) as cnt, cell_index from df group by cell_index"))
  
  #Keep first taxid per cell
  max_taxid_per_cell <- max_taxid_per_cell[!duplicated(max_taxid_per_cell$cell_index),]
  #dim(max_taxid_per_cell)
  
  #Add info to each cell
  taxid_class_per_cell <- merge(max_taxid_per_cell, df_taxid, all.x=TRUE)

  taxid_class_per_cell$cell_id <- rownames(mat)[taxid_class_per_cell$cell_index]
  taxid_class_per_cell
}




###############################################
#' Using a KRAKEN2 count matrix, produce a "kneeplot" of species
#' 
#' @param adata A Seurat object
#' @return A ggplot object
#' @export
KrakenSpeciesDistribution <- function(
    adata, 
    use_assay="kraken"
){
  strain_cnt <- adata@assays[[use_assay]]$counts
  df <- data.frame(
    taxid=rownames(strain_cnt),
    cnt=rowSums(strain_cnt)
  )
  df <- df[order(df$cnt, decreasing = TRUE),]
  df$index <- 1:nrow(df)
  ggplot(df, aes(index, cnt)) + 
    geom_line() + 
    theme_bw() +
    scale_x_log10()+
    scale_y_log10() +
    xlab("Taxonomy ID") +
    ylab("Kraken Feature count")
}



###############################################
#' Take a KRAKEN2 count matrix where the column is the taxonomyID.
#' Convert to a matrix where the columns instead are the names of each taxonomy.
#' Unused taxonomyID columns will not be kept
#' 
#' TODO: should append the name on the axis
#' 
#' @param mat description
#' @param keep_species_only description
#' @return A named count matrix
#' @export
SetTaxonomyNamesFeatures <- function(
    mat, 
    keep_species_only=TRUE
){
  
  #TODO can check if taxid_ is present for this function
  
  use_row <- as.integer(stringr::str_remove_all(colnames(mat), "taxid_"))
  #colnames(mat)
  #use_row <- which(Matrix::rowSums(mat)>0)
  
  
  taxonomizr::prepareDatabase(getAccessions=FALSE)
  taxid_class_per_cell <- as.data.frame(taxonomizr::getTaxonomy(
    use_row,
    desiredTaxa = c("phylum", "class", "order", "family", "genus","species")
  ))
  taxid_class_per_cell$taxid <- paste0("taxid_", use_row)

  if(keep_species_only) {
    taxid_class_per_cell$use_row <- use_row
    taxid_class_per_cell <- taxid_class_per_cell[!is.na(taxid_class_per_cell$species),]
    
    compressed_mat <- mat[, taxid_class_per_cell$taxid]
    colnames(compressed_mat) <- taxid_class_per_cell$species
    #compressed_mat    

  } else {
    #Find the best name to use
    taxid_class_per_cell$name <- taxid_class_per_cell$species
    torep <- is.na(taxid_class_per_cell$name)
    taxid_class_per_cell$name[torep] <- taxid_class_per_cell$genus[torep]  
    torep <- is.na(taxid_class_per_cell$name)
    taxid_class_per_cell$name[torep] <- "NA"
    
    #use_name <- paste0(use_row,"-",taxid_class_per_cell$name)
    
    compressed_mat <- mat[, taxid_class_per_cell$taxid]
    colnames(compressed_mat) <- paste0(use_row,"-",taxid_class_per_cell$name)
  }
    
  compressed_mat
}




if(FALSE){
  taxonomizr::getTaxonomy(
    0,
    desiredTaxa = c(
      "phylum", "class", "order", "family", "genus","species")
  )
  
    
}


