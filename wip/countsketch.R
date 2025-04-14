allfa <- Biostrings::readDNAStringSet("/husky/fromsequencer/240809_novaseq_wgs1/trimmed/ref10/all.fa")

map_len <- data.frame(
  longname=names(allfa),
  len=str_length(allfa)
)
map_len$name <- str_split_i(map_len$longname,"\\.", 1)
rownames(map_len) <- map_len$name



mat <- read.table("/husky/henriksson/atrandi/v2_wgs_novaseq1/countsketch_mat.csv", comment.char = "")

#mat <- read.table("/home/mahogny/github/bascet/miseqdata/countsketch_mat.csv")
mat <- read.table("/husky/henriksson/atrandi/simulated1/countsketch_mat.csv", comment.char = "")
#TODO what if we don't produce "cells" from plasmids?



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
Qnorm <- normalize_Q(Q)


if(FALSE){
  #Makes no sense!
  Qalt <- Q
  for(i in 1:nrow(Qalt)){
    Qalt[i,] <- Qalt[i,]/sum(Qalt[i,])
  }
}

mat
dim(mat)

mat[1:100,1:100]



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


pc <- prcomp(Q,
             center = FALSE,
             scale = FALSE)

pc <- prcomp(Qnorm[,celldepth>2000 & celldepth<100000],
             center = FALSE,
             scale = FALSE)  #dies

pc$x  #PC1 holds the pattern where features are really high or low


plot(
  pc$rotation[,1],
  pc$rotation[,2]
)

plot(
  pc$rotation[,3],
  pc$rotation[,4]
)


#pc$x # loadings

#pc$rotation

### screeplot; first PC might be depth related...
plot(pc$sdev)

### screeplot; ignoring first PC
plot(pc$sdev[-1])

### for each cell: scale to unit length
#cell_pos <- pc$rotation
cell_pos <- pc$rotation[,1:50]
cell_pos <- pc$rotation[,1:10]
#cell_pos <- pc$rotation[,2:10]
cell_fact <- sqrt(rowSums(cell_pos*cell_pos))
cell_norm <- cell_pos
for(i in 1:nrow(cell_pos)){
  cell_norm[i,] <- cell_norm[i,]/cell_fact[i]
}


#mean(celldepth>10000)

library(umap)
pc.umap <- umap(cell_norm)  #TODO exclude component1; only pick first 10 components?
#pc.umap <- umap(cell_norm)  #TODO exclude component1; only pick first 10 components?
pc.umap <- umap(pc$rotation[,1:10])  #TODO exclude component1; only pick first 10 components?
pc.umap <- umap(pc$rotation)  #TODO exclude component1; only pick first 10 components?

pc.umap <- umap(pc$rotation[,1:20])  #first components are more informative!

ggplot(
  data.frame(
    x=pc.umap$layout[,1], 
    y=pc.umap$layout[,2],
    ct=celltype), 
  aes(x,y,color=ct)) + geom_point()



################################################################################
############### UMAP on Q directly
################################################################################


library(umap)
pc.umap <- umap(t(Q))
#pc.umap <- umap(t(Qnorm))  #same depth for now, so can ignore issue
#pc.umap <- umap(t(Qalt))  #this makes 100% no sense

toplot <- data.frame(
  x=pc.umap$layout[,1], 
  y=pc.umap$layout[,2],
  len=map_len[celltype,]$len,
  ct=celltype)



sqldf::sqldf("select len, ct, count(*) as cnt from toplot group by ct order by len")


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
  toplot, aes(x,y,color=ct)) + 
  geom_point() 

ggplot(
  toplot[toplot$len>500000,], aes(x,y,color=len,label=ct)) + 
  geom_point() + geom_text()


table(toplot$len)
table(toplot$ct)






################################################################################
############### UMAP on Q directly -- for mock
################################################################################

head(Qnorm)
Qnorm[1:10,1:10]

mean(celldepth>100000)
mean(celldepth<30000)
keep_cells <- celldepth>200000
keep_cells <- celldepth<30000
mean(keep_cells)
sum(keep_cells)

library(umap)
#pc.umap <- umap(t(Q))
#pc.umap <- umap(t(Qnorm))  #same depth for now, so can ignore issue
#pc.umap <- umap(t(Qnorm[1:10,]))  #first features only
#pc.umap <- umap(t(Qnorm[1:20,keep_cells]), verbose=TRUE)  #first features only
pc.umap <- umap(t(norm_redQ[,keep_cells]), verbose=TRUE)  


toplot <- data.frame(
  x=pc.umap$layout[,1], 
  y=pc.umap$layout[,2],
  depth=celldepth[keep_cells]
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


reduceQ <- function(Q, numfeat){
  redQ <- matrix(data=0, ncol=ncol(Q), nrow=numfeat)
  numblock <- nrow(Q)/numfeat
  for(i in 1:numblock){
    redQ <- redQ + Q[(1:numfeat)+(numfeat*(i-1)),]*sign((i%%2)-0.5)
  }
  colnames(redQ) <- colnames(Q)
  redQ  
}
redQ <- reduceQ(Q, 1000)
norm_redQ <- normalize_Q(redQ)



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
