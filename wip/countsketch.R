allfa <- Biostrings::readDNAStringSet("/husky/fromsequencer/240809_novaseq_wgs1/trimmed/ref10/all.fa")

map_len <- data.frame(
  longname=names(allfa),
  len=str_length(allfa)
)
map_len$name <- str_split_i(map_len$longname,"\\.", 1)
rownames(map_len) <- map_len$name



#mat <- read.table("/husky/henriksson/atrandi/v2_wgs_novaseq1/countsketch_mat.csv", comment.char = "")
mat <- as.data.frame(data.table::fread("/husky/henriksson/atrandi/v2_wgs_novaseq1/countsketch_mat.csv"))

#mat <- read.table("/home/mahogny/github/bascet/miseqdata/countsketch_mat.csv")
#mat <- read.table("/husky/henriksson/atrandi/simulated1/countsketch_mat.csv", comment.char = "")
#TODO what if we don't produce "cells" from plasmids?
mat <- as.data.frame(data.table::fread("/husky/henriksson/atrandi/simulated1/countsketch_mat.csv"))
mat <- as.data.frame(data.table::fread("/husky/henriksson/atrandi/optimize_trim/countsketch_mat.csv"))


######## Break apart the input file
cellid <- mat[,1]
celldepth <- mat[,2]
Q <- t(mat[,-(1:2)])  #Each column is one cell
celltype <- str_split_i(cellid,"#",1)

colnames(Q) <- cellid
rownames(Q) <- paste0("f",1:nrow(Q))


####### Read-depth normalized version of Q. but maybe better to normalize later after decomposition?
normalize_Q <- function(Q){
  Qnorm <- Q
  for(i in 1:ncol(Q)){
    Qnorm[,i] <- Qnorm[,i]/(10000*celldepth[i]) #factor to avoid too small numbers
  }
  Qnorm
}


var_normalize_Q <- function(Q){
  Qnorm <- Q
  for(i in 1:ncol(Q)){
#    print(sqrt(Qnorm[,i]))
    Qnorm[,i] <- Qnorm[,i]/max(abs(Qnorm[,i])) #factor to avoid too small numbers
#    Qnorm[,i] <- Qnorm[,i]/(var(Qnorm[,i])) #factor to avoid too small numbers
#    Qnorm[,i] <- Qnorm[,i]/sqrt(var(Qnorm[,i])) #factor to avoid too small numbers
  }
  Qnorm
}


reduceQ <- function(Q, numfeat){
  redQ <- matrix(data=0, ncol=ncol(Q), nrow=numfeat)
  numblock <- nrow(Q)/numfeat
  for(i in 1:numblock){
    redQ <- redQ + Q[(1:numfeat)+(numfeat*(i-1)),]*sign((i%%2)-0.5)
  }
  colnames(redQ) <- colnames(Q)
  redQ  
}

Qnorm <- normalize_Q(Q)
redQ <- reduceQ(Q, 500) 
norm_redQ <- normalize_Q(redQ) #100 works for mock, if normalized

var_norm_redQ <- var_normalize_Q(redQ) #100 works for mock, if normalized


### TODO compare to KRAKEN!
get_avg_prof <- function(Q){
  all_list <- list()
  for(ct in unique(celltype)) {
    all_list[[ct]] <- data.frame(
      ct=ct,
      feat=1:nrow(Q),
      lev=rowMeans(Q[,celltype==ct,drop=FALSE])  
    )
  }
  do.call(rbind, all_list)
}
avg_prof <- get_avg_prof(redQ)
ggplot(avg_prof, aes(feat, lev, color=ct)) + geom_line()




#https://www.rdocumentation.org/packages/aroma.light/versions/3.2.0/topics/wpca


colMeans(Q)

plot(sort(celldepth))
plot(sort(rowMeans(Q))) #around 0 ... but few features are "random" near 0. expected??
plot(sort(colMeans(Q))) #per cell: around 0

plot(colMeans(Q), celldepth)


#one feature over all cells
plot(sort(Q[1,]))
plot(sort(Q[2,]))
plot(Q[1,], celldepth)
plot(Q[2,], celldepth)

#one cell, all features
plot(sort(Q[,1]))
plot(sort(Q[,2]))



plot(Qnorm[1,], Qnorm[2,])
# 
# 
# pc <- prcomp(Q,
#              center = FALSE,
#              scale = FALSE)
# 
# pc <- prcomp(Qnorm[,celldepth>2000 & celldepth<100000],
#              center = FALSE,
#              scale = FALSE)  #dies
# 
# pc$x  #PC1 holds the pattern where features are really high or low
# 
# 
# plot(
#   pc$rotation[,1],
#   pc$rotation[,2]
# )
# 
# plot(
#   pc$rotation[,3],
#   pc$rotation[,4]
# )
# 
# 
# #pc$x # loadings
# 
# #pc$rotation
# 
# ### screeplot; first PC might be depth related...
# plot(pc$sdev)
# 
# ### screeplot; ignoring first PC
# plot(pc$sdev[-1])
# 
# ### for each cell: scale to unit length
# #cell_pos <- pc$rotation
# cell_pos <- pc$rotation[,1:50]
# cell_pos <- pc$rotation[,1:10]
# #cell_pos <- pc$rotation[,2:10]
# cell_fact <- sqrt(rowSums(cell_pos*cell_pos))
# cell_norm <- cell_pos
# for(i in 1:nrow(cell_pos)){
#   cell_norm[i,] <- cell_norm[i,]/cell_fact[i]
# }


#mean(celldepth>10000)
# 
# library(umap)
# pc.umap <- umap(cell_norm)  #TODO exclude component1; only pick first 10 components?
# #pc.umap <- umap(cell_norm)  #TODO exclude component1; only pick first 10 components?
# pc.umap <- umap(pc$rotation[,1:10])  #TODO exclude component1; only pick first 10 components?
# pc.umap <- umap(pc$rotation)  #TODO exclude component1; only pick first 10 components?
# 
# pc.umap <- umap(pc$rotation[,1:20])  #first components are more informative!
# 
# ggplot(
#   data.frame(
#     x=pc.umap$layout[,1], 
#     y=pc.umap$layout[,2],
#     ct=celltype), 
#   aes(x,y,color=ct)) + geom_point()



################################################################################
############### UMAP on Q directly - general
################################################################################

#order metadata according to input data
metadata <- readRDS("/husky/henriksson/atrandi/v2_wgs_novaseq1/metadata.RDS")
metadata$cellid <- stringr::str_remove(rownames(metadata),"BASCET_")
rownames(metadata) <- metadata$cellid
metadata <- metadata[cellid,]

keep_cells <- 1:1000

keep_cells <- celldepth>40000 & celldepth<50000 & metadata$genus!="Cereibacter" & metadata$genus!="Homo"
sum(keep_cells)

library(umap)
#pc.umap <- umap(t(Q), verbose=TRUE)
pc.umap <- umap(t(Qnorm), verbose=TRUE)  #same depth for now, so can ignore issue
pc.umap <- umap(t(norm_redQ[,keep_cells]), verbose=TRUE)  #same depth for now, so can ignore issue
pc.umap <- umap(t(var_norm_redQ[,keep_cells]), verbose=TRUE)  #same depth for now, so can ignore issue


toplot <- data.frame(
  cellid=cellid[keep_cells],
  x=pc.umap$layout[,1], 
  y=pc.umap$layout[,2],
  depth=celldepth[keep_cells],
  var=sqrt(apply(norm_redQ[,keep_cells],2,var))
)

ggplot(
  merge(toplot, metadata), aes(x,y,color=genus)) + 
  geom_point() 

ggplot(
  toplot, aes(x,y,color=log10(depth))) + 
  geom_point() 

ggplot(
  toplot, aes(x,y,color=log(var))) + 
  geom_point() 



egg::ggarrange(
  ggplot(
    toplot, aes(x,y,color=len,label=ct)) + 
    geom_point() + geom_text(),
  ggplot(
    toplot[toplot$len>500000,], aes(x,y,color=len,label=ct)) + 
    geom_point() + geom_text()
)

## just color
egg::ggarrange(
  ggplot(
    toplot, aes(x,y,color=ct)) + 
    geom_point(), 
  ggplot(
    toplot[toplot$len>2000000,], aes(x,y,color=ct)) + 
    geom_point() 
)



ggplot(
  toplot[toplot$len>500000,], aes(x,y,color=len,label=ct)) + 
  geom_point() + geom_text()


#table(toplot$len)
#table(toplot$ct)






################################################################################
############### UMAP on Q directly -- for mock
################################################################################

head(Qnorm)
Qnorm[1:10,1:10]

mean(celldepth>100000)
mean(celldepth<30000)
keep_cells <- celldepth>200000
keep_cells <- celldepth<25000
mean(keep_cells)
sum(keep_cells)

library(umap)
#pc.umap <- umap(t(Q))
#pc.umap <- umap(t(Qnorm))  #same depth for now, so can ignore issue
#pc.umap <- umap(t(Qnorm[1:10,]))  #first features only
#pc.umap <- umap(t(Qnorm[1:20,keep_cells]), verbose=TRUE)  #first features only
#pc.umap <- umap(t(norm_redQ[,keep_cells]), verbose=TRUE)  
#pc.umap <- umap(t(redQ[,keep_cells]), verbose=TRUE)
#pc.umap <- umap(t(Q[1:20,keep_cells]), verbose=TRUE)
pc.umap <- umap(t(redQ), verbose=TRUE)
pc.umap <- umap(t(norm_redQ), verbose=TRUE)  #normalization helps a lot if only 100 features



toplot <- data.frame(
  x=pc.umap$layout[,1], 
  y=pc.umap$layout[,2],
  len=map_len[celltype,]$len,
  depth=celldepth,
  ct=celltype
)

egg::ggarrange(
  ggplot(
    toplot, aes(x,y,color=depth)) + 
    geom_point() ,
  ggplot(
    toplot, aes(x,y,color=ct)) + 
    geom_point()
)


ggplot(
  toplot, aes(x,y,color=log10(depth))) + 
  geom_point() 


#machine precision issues?
#solve umap with division on the fly?










########### Hadamard-like matrix

#n

if(FALSE){
  numfeat <- 10
  hadmat <- matrix(data = 0, ncol=nrow(Q), nrow=numfeat)
  for(i in 1:numfeat){
    hadmat[i,]
  }
}






################################################################################
############### Seurat version
################################################################################


#put norm_redQ as PC

#RunPCA
#https://satijalab.org/seurat/reference/runpca

#https://github.com/satijalab/seurat/blob/HEAD/R/dimensional_reduction.R

## NN graph takes a lot of time for regular umap. seurat better?



eduction.data <- CreateDimReducObject(
  embeddings = cell.embeddings,
  loadings = feature.loadings,
  assay = assay,
  stdev = sdev,
  key = reduction.key,
  misc = list(total.variance = total.variance)
)







pbmc <- CreateSeuratObject(counts = norm_redQ[,keep_cells], project = "pbmc3k", min.cells = 0, min.features = 0)
#Warning: Data is of class matrix. Coercing to dgCMatrix.

pbmc <- FindNeighbors(pbmc, dims = 1:10, reduction = )
pbmc <- FindClusters(pbmc, resolution = 0.5)


pbmc <- RunUMAP(pbmc, dims = 1:10)



#TODO try cells with low coverage!!





#FindNeighbors
#https://github.com/satijalab/seurat/blob/95ed7691a25365338848e3532733454915b475da/R/clustering.R#L652










RunQ.default <- function(
    object,
    assay = NULL,
    reduction.key = "LSI_",
    verbose = TRUE,
    ...
) {
  if (is.null(x = rownames(x = object))) {
    rownames(x = object) <- seq_len(length.out = nrow(x = object))
  }
  if (is.null(x = colnames(x = object))) {
    colnames(x = object) <- seq_len(length.out = ncol(x = object))
  }
  
  
  
  n <- min(n, (ncol(x = object) - 1))
  if (verbose) {
    message("Running Q")
  }
  components <- irlba(A = t(x = object), nv = n, work = irlba.work, tol = tol)
  feature.loadings <- components$v
  sdev <- components$d / sqrt(x = max(1, nrow(x = object) - 1))
  cell.embeddings <- components$u
  if (scale.embeddings) {
    if (verbose) {
      message("Scaling cell embeddings")
    }
    embed.mean <- apply(X = cell.embeddings, MARGIN = 2, FUN = mean)
    embed.sd <- apply(X = cell.embeddings, MARGIN = 2, FUN = sd)
    norm.embeddings <- t((t(cell.embeddings) - embed.mean) / embed.sd)
    if (!is.null(x = scale.max)) {
      norm.embeddings[norm.embeddings > scale.max] <- scale.max
      norm.embeddings[norm.embeddings < -scale.max] <- -scale.max
    }
  } else {
    norm.embeddings <- cell.embeddings
  }
  rownames(x = feature.loadings) <- rownames(x = object)
  colnames(x = feature.loadings) <- paste0(
    reduction.key, seq_len(length.out = n)
  )
  rownames(x = norm.embeddings) <- colnames(x = object)
  colnames(x = norm.embeddings) <- paste0(
    reduction.key, seq_len(length.out = n)
  )
  reduction.data <- CreateDimReducObject(
    embeddings = norm.embeddings,
    #loadings = feature.loadings,
    assay = assay,
    stdev = sdev,
    key = reduction.key,
    misc = components
  )
  return(reduction.data)
}


RunQ.Assay <- function(
    object,
    assay = NULL,
    features = NULL,
    reduction.key = "LSI_",
    verbose = TRUE,
    ...
) {
  features <- SetIfNull(x = features, y = VariableFeatures(object = object))
  data.use <- GetAssayData(
    object = object,
    layer = "data"
  )[features, ]
  reduction.data <- RunQ(
    object = data.use,
    assay = assay,
    features = features,
    n = n,
    reduction.key = reduction.key,
    scale.max = scale.max,
    verbose = verbose,
    ...
  )
  return(reduction.data)
}


RunQ.StdAssay <- function(
    object,
    assay = NULL,
    features = NULL,
    reduction.key = "LSI_",
    verbose = TRUE,
    ...
) {
  RunQ.Assay(
    object = object,
    assay = assay,
    features = features,
    reduction.key = reduction.key,
    verbose = verbose,
    ...
  )
}

RunQ.Seurat <- function(
    object,
    assay = NULL,
    features = NULL,
    reduction.key = "LSI_",
    reduction.name = "lsi",
    verbose = TRUE,
    ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  assay.data <- object[[assay]]
  reduction.data <- RunQ(
    object = assay.data,
    assay = assay,
    features = features,
    reduction.key = reduction.key,
    verbose = verbose,
    ...
  )
  object[[reduction.name]] <- reduction.data
  return(object)
}

