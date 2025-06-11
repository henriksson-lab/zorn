library(stringr)
library(Biostrings)


################################################################################ 
###### Make a mock community, with plasmids, roughly uniform distribution ######
################################################################################


insert_max_size <- 800
num_cells <- 50
num_shards <- 20

fa_rootdir <- "/husky/fromsequencer/240809_novaseq_wgs1/trimmed/ref10/separate"
list_genomes <- unique(str_split_i(list.files(fa_rootdir)," ",1))


pos_index <- 1
for(cur_shard in 1:num_shards) {
  allreads_df <- list()
  print(paste("shard", cur_shard))
  for(cur_genome in list_genomes) {
    
    list_fa <- list.files(fa_rootdir)
    list_fa <- list_fa[stringr::str_starts(list_fa,cur_genome)]
    print(list_fa)
    
    allfa <- Biostrings::readDNAStringSet(file.path(fa_rootdir, list_fa)) 
    allfa_len <- str_length(allfa)
    allfa_str <- as.character(allfa)
    
    #### For each cell
    for(cell_i in 1:num_cells) {
      
      #num_reads <- 1000  #Randomize on log scale
      num_reads <- round(10**runif(1,3,4.5))
      
      #Decide chromosome to use
      use_chrom <- sample(1:length(allfa), prob = allfa_len, replace = TRUE, size = num_reads)
      
      
      
      readlen <- 150
      qread <- str_flatten(rep("F",readlen))
      
      ### For each chromosome to include
      #tab_use_chrom <- table(use_chrom)
      for(cur_chrom in unique(use_chrom)){
        get_num_reads <- sum(use_chrom==cur_chrom)
        
        ### Generate reads from this chromosome
        ch <- allfa_str[cur_chrom]
        totlen <- str_length(ch)
        start_pos <- round(runif(get_num_reads, min=1, max=totlen-insert_max_size))
        end_pos <- start_pos + round(runif(get_num_reads, min=400, max=insert_max_size))
        
        all_inserts <- str_sub(ch, start_pos, end_pos)
        
        insert_r1 <- str_sub(all_inserts, 1, readlen)
        insert_r2_rc <- str_sub(all_inserts, str_length(all_inserts)-(readlen-1), str_length(all_inserts))
        insert_r2 <- as.character(reverseComplement(DNAStringSet(insert_r2_rc)))
        
        df <- data.frame(
          name=sprintf("cell%05d",pos_index),
          pos1=1,
          pos2=1,
          r1=insert_r1,
          r2=insert_r2,
          q1=qread,
          q2=qread,
          umi=""
        )  
        
        allreads_df[[pos_index]] <- df
        pos_index <- pos_index + 1
        
      }
    }
  }
  
  ## Flatten list, save
  allreads_flat <- do.call(rbind, allreads_df)  
  
  fname <- paste0("/husky/henriksson/atrandi/simulated4/filtered.",cur_shard,".tirp")  
  
  print("writing fastq")
  write.table(allreads_flat, fname, row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)   
  
  print("compressing")
  system(paste0("bgzip --threads 15 -f ",fname))
  print("indexing")
  system(paste0("tabix -f -p bed ",fname, ".gz"))
  
}












