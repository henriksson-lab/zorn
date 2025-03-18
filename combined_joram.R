
list_adata <- list()

list_index <- c("10min","5min")
for(curlib in 1:2){
  print(curlib)
  bascetRoot <- file.path("/husky/henriksson/atrandi/rnaseq3",curlib)
  adata_f <- FragmentsToSignac(file.path(bascetRoot,"fragments.1.tsv.gz"))  
  adata <- CreateSeuratObject(CountGrangeFeatures(adata_f, grange_gene))
  adata$index <- list_index[curlib]
  
  ##divide library
  adata$well_r1 <- stringr::str_split_i(colnames(adata),"_",2)
  AtoH <- c("A","B","C","D","E","F","G","H")
  libname <- c(rep("s1",8),rep("s3",8),rep("s2",8))  #note order
  names(libname) <- c(paste0(AtoH,1), paste0(AtoH,2), paste0(AtoH,3))
  adata$libname <- unname(libname[adata$well_r1])
  
  
  #s 1: 0.2 U/µL papA, 1-minute incubation
  #s 2: 0.2 U/µL papA, 2-minute incubation
  #s 3: 0.2 U/µL papA, 1-minute incubation   no dnase treatment

  list_adata[[curlib]] <- adata  
}

adata <- merge(list_adata[[1]], y = list_adata[[2]], add.cell.ids = c("a", "b"), project = "PBMC12K")
adata$index_lib <- paste(adata$index, adata$libname)

#################### knee plots



### for each library
all_df <- NULL
for(curlib in unique(adata$index_lib)){
  df <- data.frame(cnt=sort(adata$nCount_RNA[adata$index_lib==curlib], decreasing = TRUE))
  df$rank <- 1:nrow(df)
  df$lib=curlib
  all_df <- rbind(all_df, df)
}
ggplot(all_df,aes(rank, cnt, color=lib, group=lib)) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() + theme_bw()
ggsave("~/kemal/kneeplot.pdf", width = 6, height = 6)

################################################


adata <- adata[,adata$nCount_RNA>10] #meaningless below for certain!!

adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(adata), 10)
top10
#[1] "npr"              "glmM"             "rpoE"             "tsf"              "rplB"             "hupB"            
#[7] "bamA"             "lpxC"             "ssrA"             "rna-BZ22-RS22510"

all.genes <- rownames(adata)
adata <- ScaleData(adata, features = all.genes)
adata <- RunPCA(adata, features = VariableFeatures(object = adata))
adata <- FindNeighbors(adata, dims = 1:10)
adata <- RunUMAP(adata, dims = 1:10)

adata <- FindClusters(adata, resolution = 0.2)

DimPlot(adata, reduction = "umap")

DimPlot(adata, reduction = "umap", group.by = "index")
ggsave("~/kemal/umap_by_time.pdf", width = 6, height = 6)

DimPlot(adata, reduction = "umap", group.by = "libname")
ggsave("~/kemal/umap_by_libname.pdf", width = 6, height = 6)


FeaturePlot(adata, features=top10[1:8])

FeaturePlot(adata, features=c("ssrA","rna-BZ22-RS22510","BZ22-RS01385","fghA"))
ggsave("~/kemal/umap_example_genes.pdf", width = 12, height = 12)

#FeaturePlot(adata, features=c("npr"))

################
################
################




alldf <- NULL
for(curlib in unique(adata$index_lib)){
  df <- as.data.frame(merge(
    data.frame(
      Name=rownames(adata),
      cnt=rowSums(adata[,adata$index_lib==curlib]@assays$RNA@counts),
      lib=curlib
    ),
    as.data.frame(grange_gene)[,c("Name","gene_biotype")]
  ))
  #as.data.frame(df)
  alldf <- rbind(
    alldf,
    sqldf::sqldf("select gene_biotype, sum(cnt) as cnt, lib from df group by gene_biotype")
  )
  #protein_coding is the relevant one  
}
sum_cnt_type <- reshape2::acast(data = alldf, gene_biotype~lib, value.var = "cnt")

sum_cnt_type_norm <- sum_cnt_type
for(i in 1:ncol(sum_cnt_type_norm)){
  sum_cnt_type_norm[,i] <- round(digits = 3,sum_cnt_type_norm[,i]/sum(sum_cnt_type_norm[,i]))
}

sum_cnt_type
sum_cnt_type_norm  ### I do not believe this data




