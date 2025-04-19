#note scann seems better than annoy
#https://ann-benchmarks.com/glove-100-angular_10_angular.html

allfa <- Biostrings::readDNAStringSet("/husky/fromsequencer/240809_novaseq_wgs1/trimmed/ref10/all.fa")

map_len <- data.frame(
  longname=names(allfa),
  len=str_length(allfa)
)
map_len$name <- str_split_i(map_len$longname,"\\.", 1)
rownames(map_len) <- map_len$name


#mat <- read.table("/husky/henriksson/atrandi/v2_wgs_novaseq1/countsketch_mat.csv", comment.char = "")
mat <- as.data.frame(data.table::fread("/husky/henriksson/atrandi/v2_wgs_novaseq1/countsketch_mat.csv"))
mat <- as.data.frame(data.table::fread("/husky/henriksson/atrandi/v3_wgs_novaseq1/countsketch_mat.csv"))
mat <- as.data.frame(data.table::fread("/husky/henriksson/atrandi/v2_wgs_saliva1/countsketch_mat.csv"))
mat <- as.data.frame(data.table::fread("/husky/henriksson/atrandi/v4_wgs_saliva1/countsketch_mat.csv"))


#mat <- read.table("/home/mahogny/github/bascet/miseqdata/countsketch_mat.csv")
#mat <- read.table("/husky/henriksson/atrandi/simulated1/countsketch_mat.csv", comment.char = "")
#TODO what if we don't produce "cells" from plasmids?
mat <- as.data.frame(data.table::fread("/husky/henriksson/atrandi/simulated1/countsketch_mat.csv"))
mat <- as.data.frame(data.table::fread("/husky/henriksson/atrandi/optimize_trim/countsketch_mat.csv"))  
#262M reads for 240903_wgs_atcc2_miseq/raw ?? 5422604 lines in asfq. or 1.3M reads
#dim(mat)


######## Break apart the input file
cellid <- mat[,1]
celldepth <- mat[,2]
Q <- t(mat[,-(1:2)])  #Each column is one cell
celltype <- str_split_i(cellid,"#",1)

colnames(Q) <- cellid
rownames(Q) <- paste0("f",1:nrow(Q))

####### QC
sum(celldepth) #this is the number of kmers... 
sum(celldepth)/(150+150-40) #number of reads, approx


####### Read-depth normalized version of Q. using kmer count
normalize_Q <- function(Q){
  Qnorm <- Q #normalize column-wise
  for(i in 1:ncol(Q)){
    Qnorm[,i] <- Qnorm[,i]*(10000/celldepth[i]) #factor to avoid too small numbers   ... beep! change!!!
  }
  Qnorm
}

####### Read-depth normalized version of Q. size factor
sf_normalize_Q <- function(Q){
  Qnorm <- Q #normalize column-wise
  for(i in 1:ncol(Q)){
    Qnorm[,i] <- Qnorm[,i]*(10000/sum(Qnorm[,i])) #factor to avoid too small numbers
  }
  Qnorm
}

####### Read-depth normalized version of Q. normalize to unit length
l2_normalize_Q <- function(Q){
  Qnorm <- Q #normalize column-wise
  for(i in 1:ncol(Q)){
    Qnorm[,i] <- Qnorm[,i]/sqrt(sum(Qnorm[,i]*Qnorm[,i])) #factor to avoid too small numbers
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


remove_comp1 <- function(Q){
  #Regress out first component.
  #in the order of 0! (this is good) so meaningless
  m <- rowMeans(Q)
  len_m <- sqrt(sum(m*m))
  for(i in 1:ncol(Q)){
    print(sum(Q[,i]*m)/len_m)
    Q[,i] <- Q[,i] - sum(Q[,i]*m)/len_m*m
  }
  Q
}
#norm_redQ_no1 <- remove_comp1(norm_redQ)
#norm_redQ_no1[1:5,1:5]
#remove_comp1(remove_comp1(norm_redQ))[1:5,1:5]
#norm_redQ[1:5,1:5]

#var_norm_redQ <- var_normalize_Q(redQ) #100 works for mock, if normalized


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
#avg_prof <- get_avg_prof(redQ)
#ggplot(avg_prof, aes(feat, lev, color=ct)) + geom_line()




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



################################################################################
############### UMAP on Q directly - general with kraken comparison
################################################################################

#order metadata according to input data
metadata <- readRDS("/husky/henriksson/atrandi/v2_wgs_novaseq1/metadata.RDS")
metadata$cellid <- stringr::str_remove(rownames(metadata),"BASCET_")
rownames(metadata) <- metadata$cellid
metadata <- metadata[cellid,]

keep_cells <- 1:1000

keep_cells <- celldepth>40000 & celldepth<50000 & metadata$genus!="Cereibacter" & metadata$genus!="Homo"
keep_cells <- celldepth>200000
keep_cells <- celldepth>2000000
sum(keep_cells)

library(umap)
#pc.umap <- umap(t(Q), verbose=TRUE)
pc.umap <- umap(t(Qnorm), verbose=TRUE)  #same depth for now, so can ignore issue
pc.umap <- umap(t(Qnorm[,keep_cells]), verbose=TRUE)  #same depth for now, so can ignore issue
pc.umap <- umap(t(norm_redQ[,keep_cells]), verbose=TRUE)  #same depth for now, so can ignore issue
pc.umap <- umap(t(sign(norm_redQ[,keep_cells])), verbose=TRUE)  #same depth for now, so can ignore issue
pc.umap <- umap(t(sign(Q[,keep_cells])), verbose=TRUE)  #same depth for now, so can ignore issue
#pc.umap <- umap(t(var_norm_redQ[,keep_cells]), verbose=TRUE)  #same depth for now, so can ignore issue


toplot <- data.frame( ############ how is xantham affecting?
  cellid=cellid[keep_cells],
  x=pc.umap$layout[,1], 
  y=pc.umap$layout[,2],
  depth=celldepth[keep_cells]
)

ggplot(
  toplot, aes(x,y,color=log10(depth))) + 
  geom_point() 


ggplot(
  merge(toplot, metadata), aes(x,y,color=genus)) + 
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


library(plotly)
toplot2 <- merge(toplot, metadata)
fig <- plot_ly(type = 'scatter', mode = 'markers') 
fig <- fig %>%
  add_trace(
    x = toplot2$x, 
    y = toplot2$y,
    text = toplot2$genus,
    hoverinfo = 'text',
    marker = list(color='green'),
    showlegend = F
  )
fig



library(plotly)
toplot2 <- merge(toplot, metadata)
fig <- plot_ly(toplot2, x=~x, y=~y, split=~genus, type = 'scatter', mode = 'markers') 
fig

library(plotly)
toplot2 <- merge(toplot, metadata)
fig <- plot_ly(toplot2[toplot2$genus!="Homo",], x=~x, y=~y, split=~genus, type = 'scatter', mode = 'markers') 
fig



################################################################################
############### UMAP on Q directly -- for simulated
################################################################################
library(umap)

pc.umap <- umap(t(sign(norm_redQ)), verbose=TRUE)  


toplot <- data.frame(
  x=pc.umap$layout[,1], 
  y=pc.umap$layout[,2],
  len=map_len[celltype,]$len,
  depth=celldepth,
  ct=celltype
)

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

################################################################################
############### UMAP on Q directly -- for mock
################################################################################

#note: depth here is in number of kmers i.e. bases

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
#pc.umap <- umap(t(Qnorm[1:20,keep_cells]), verbose=TRUE)  #first features only
pc.umap <- umap(t(norm_redQ[,keep_cells]), verbose=TRUE)  
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

toplot <- data.frame(
  x=pc.umap$layout[,1], 
  y=pc.umap$layout[,2],
  len=map_len[celltype,]$len[keep_cells],
  depth=celldepth[keep_cells],
  ct=celltype[keep_cells]
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



CreateSeuratObjectWithReduction(Q)



#Warning: Data is of class matrix. Coercing to dgCMatrix.
#pbmc <- FindNeighbors(pbmc, dims = 1:10, reduction = "kmersketch", annoy.metric = "cosine")
#
#pbmc <- FindClusters(pbmc, resolution = 0.5)
#TODO try cells with low coverage!!
#FindNeighbors
#https://github.com/satijalab/seurat/blob/95ed7691a25365338848e3532733454915b475da/R/clustering.R#L652

library(future)
plan("multicore", workers = 10)

keep_cells <- celldepth>2000000
keep_cells <- celldepth>200000
keep_cells <- celldepth>0
sum(keep_cells)

pbmc <- CreateSeuratObject(counts = Q[1:2,keep_cells,drop=FALSE], project = "pbmc3k", min.cells = 0, min.features = 0)

if(FALSE){
  pbmc$ct <- str_split_i(colnames(Q),"#",1)   #celltype
}

if(FALSE){
  #Add kraken metadata
  metadata <- readRDS("/husky/henriksson/atrandi/v2_wgs_novaseq1/metadata.RDS")
  metadata$cellid <- stringr::str_remove(rownames(metadata),"BASCET_")
  rownames(metadata) <- metadata$cellid
  metadata <- metadata[cellid,]
  pbmc$ct <- metadata[colnames(pbmc),]$genus
  pbmc$depth <- celldepth
}



pbmc@reductions[["kmersketch"]] <- CreateDimReducObject(
#  embeddings = t(norm_redQ[,keep_cells]),
#  embeddings = t(norm_redQ[,keep_cells]),
#  embeddings = t(Q[,keep_cells]),
  embeddings = sign(t(Q[,keep_cells])), 
#  embeddings = sign(t(redQ[,keep_cells])), 
  key = "kmersketch",
  assay = "RNA"
)
pbmc <- RunUMAP(pbmc, dims = 1:ncol(pbmc@reductions$kmersketch@cell.embeddings), reduction = "kmersketch")  #Searching Annoy index using 1 thread, search_k = 3000 ; can do more
DimPlot(pbmc, group.by = "ct")
DimPlot(pbmc[,pbmc$ct!="Homo"], group.by = "ct")
DimPlot(pbmc)
FeaturePlot(pbmc, features = "depth")


##### not convinced that this does the right thing
pbmc <- FindNeighbors(pbmc, reduction = "kmersketch", annoy.metric = "cosine")
pbmc <- FindClusters(pbmc, resolution = 0.05)#, reduction = "kmersketch")
DimPlot(pbmc)#, group.by = "ct")



################################################################################
############### PCA analysis
################################################################################



#fit <- princomp(t(Qnorm))#[,1:100])
#screeplot(fit, npcs = 100)
#screeplot(fit, npcs = 24, type = "lines")

#fit <- princomp(norm_redQ)#[,1:100])
#screeplot(fit, npcs = 100)


pca <- prcomp(t(Qnorm), scale = FALSE)
screeplot(pca, npcs = 10)

pca <- prcomp(t(norm_redQ), scale = FALSE)
screeplot(pca, npcs = 10)


################################################################################
############### ICA analysis
################################################################################


library(ica)

imod <- icafast(t(norm_redQ), 20)

imod$M
imod$S

plot(imod$S[,1],imod$S[,2])


acy(Bmat, imod$M)
cor(Amat, imod$S)

#, center = TRUE, maxit = 100, tol = 1e-6, Rmat = diag(nc),
#        alg = "par", fun = "logcosh", alpha = 1)


pca <- prcomp(t(norm_redQ), scale = FALSE)
screeplot(pca, npcs = 10)


pc.umap <- umap(imod$S, verbose=TRUE)  #normalization helps a lot if only 100 features
pc.umap


toplot <- data.frame(
  x=pc.umap$layout[,1], 
  y=pc.umap$layout[,2],
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



toplot <- data.frame(
  x=pc.umap$layout[,1], 
  y=pc.umap$layout[,2],
  depth=celldepth,
  ct=metadata[cellid,]$species #celltype
)


ggplot(
  toplot[sample(1:nrow(toplot),size = 10000),], aes(x,y,color=ct)) + 
  geom_point()







################################################################################
############### Compare to known genomes
################################################################################


mat <- as.data.frame(data.table::fread("/husky/henriksson/atrandi/simulated1/allgenome/countsketch_mat.csv"))  #each genome
W <- t(mat[,-(1:2)])  #Each column is one cell
colnames(W) <- str_split_i(mat[,1], stringr::fixed("."), 1)
rownames(W) <- paste0("f",1:nrow(W))


proj_QW <- t(W) %*% Q
proj_QW <- t(W) %*% Qnorm


#norm_proj_QW <- normalize_Q(proj_QW)
norm_proj_QW <- sf_normalize_Q(proj_QW)
norm_proj_QW <- l2_normalize_Q(proj_QW)


#normalize_Q
#proj_QW[1:10,1:10]

#### projected vectors are NOT sf-normalized!
plot(sort(as.double(colSums(proj_QW))))

#### projected vectors are NOT l2-normalized!
plot(sort(sqrt(as.double(colSums(proj_QW*proj_QW)))))


dim(W)
dim(Q)

keep_cells <- 1:1000

keep_cells <- celldepth>100000
sum(keep_cells)

pc.umap <- umap(t(proj_QW[,keep_cells]), verbose=TRUE)  
pc.umap <- umap(t(norm_proj_QW[,keep_cells]), verbose=TRUE)  


toplot <- data.frame(
  x=pc.umap$layout[,1], 
  y=pc.umap$layout[,2],
  depth=celldepth[keep_cells],
  ct=metadata[cellid,]$species[keep_cells]
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




toplot <- data.frame(
  x=pc.umap$layout[,1], 
  y=pc.umap$layout[,2],
  p=proj_QW[4,keep_cells]
)
ggplot(toplot, aes(x,y,color=p)) + geom_point() 
rownames(proj_QW)
#NZ_LR881938 in tail
#CP006777 away from tail





################################################################################
############### SVD??
################################################################################

redq <- Qnorm[,keep_cells]
dim(redq)

#https://rdrr.io/cran/irlba/man/irlba.html

#


usvd <- irlba::irlba(redq,nv=50)
#??irlba

princ.comps <- usvd$u %*% diag(usvd$d)
#cell.embeddings <- usvd$u %*% diag(usvd$d)
plot(log10(usvd$d))

#https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
#If 𝐗=𝐔𝐒𝐕⊤, then the columns of 𝐕 are principal directions/axes (eigenvectors)
#Columns of 𝐔𝐒 are principal components ("scores").


#TODO: compare V with pure genome vectors. similarity?


l2norm_svd <- t(l2_normalize_Q(t(usvd$v)))
l2norm_svd <- t(l2_normalize_Q(t(usvd$v[,-1])))

pc.umap <- umap(usvd$v, verbose=TRUE)  
pc.umap <- umap(usvd$v[,-1], verbose=TRUE)  
pc.umap <- umap(l2norm_svd, verbose=TRUE)  


plot(sort(rowSums(usvd$v*usvd$v))) #l2 
plot(sort(rowSums(l2norm_svd*l2norm_svd)))


toplot <- data.frame(
  x=pc.umap$layout[,1], 
  y=pc.umap$layout[,2],
  depth=celldepth[keep_cells],
  ct=metadata[cellid,]$species[keep_cells]
)


egg::ggarrange(
  ggplot(
    toplot, aes(x,y,color=depth)) + 
    geom_point() ,
  ggplot(
    toplot, aes(x,y,color=ct)) + 
    geom_point()
)


#https://github.com/bnprks/BPCells


if(FALSE){
  #this library is for large scale single-cell analysis. has svd etc
  remotes::install_github("bnprks/BPCells/r")
  library(BPCells)
  svd.function <- function(A, nv, ...) BPCells::svds(A=A, k = nv)
  #https://github.com/bnprks/BPCells/blob/b351734c9a3122f799a2e6ac5b82b77b75d05d5d/r/src/bpcells-cpp/matrixIterators/SVD.cpp#L46
  #support for on-disk SVD etc
}





npcs <- 10
pca.results <- svd.function(A = t(x = object), nv = npcs)
feature.loadings <- pca.results$v
sdev <- pca.results$d/sqrt(max(1, ncol(object) - 1))
if (weight.by.var) {
  cell.embeddings <- pca.results$u %*% diag(pca.results$d)
} else {
  cell.embeddings <- pca.results$u
}


