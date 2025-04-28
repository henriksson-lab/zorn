library(stringr)
library(Biostrings)

############### Load fasta to simulate
allfa <- Biostrings::readDNAStringSet("/husky/fromsequencer/240809_novaseq_wgs1/trimmed/ref10/all.fa")


readlen <- 150
qread <- str_flatten(rep("F",readlen))

simulate_one_cell <- function(ch, genome_name, num_read=1000, insert_max_size=800) {
  totlen <- str_length(ch)
  start_pos <- round(runif(num_read, min=1, max=totlen-insert_max_size))
  end_pos <- start_pos + round(runif(num_read, min=400, max=insert_max_size))
  
  all_inserts <- str_sub(ch, start_pos, end_pos)
  
  insert_r1 <- str_sub(all_inserts, 1, readlen)
  insert_r2_rc <- str_sub(all_inserts, str_length(all_inserts)-(readlen-1), str_length(all_inserts))
  insert_r2 <- as.character(reverseComplement(DNAStringSet(insert_r2_rc)))
  
  data.frame(
    name=genome_name, 
    pos1=1,
    pos2=1,
    r1=insert_r1,
    r2=insert_r2,
    q1=qread,
    q2=qread,
    umi=""
  )  
}


finalize_file <- function(outf){
  system(paste0("bgzip -f ",outf))
  system(paste0("tabix -f -p bed ",outf, ".gz"))  
}

################################################################################ 
############# Make a mock community that is similar to measured mock ###########
################################################################################

#Load empirical distribution of genome depths
df_align_bc_dist <- readRDS("/husky/henriksson/atrandi/v4_wgs_novaseq3/align_metadata.RDS")

#https://www.atcc.org/products/msa-2003 each strain name starts with different word

allfa_main <- allfa[str_length(allfa)>500000] #skip plasmids
names(allfa_main)
names(allfa_main) <- stringr::str_split_fixed(names(allfa_main)," ",3)[,2]
unique(names(allfa_main))

# * use count distribution from mock3
# * use kraken inferred species
# * 20 random output TIRP files

num_shard <- 20
insert_max_size <- 800



print("Cleaning")
for(i in 1:num_shard){
  print(i)
  fname <- paste0("/husky/henriksson/atrandi/simulated2/filtered.",i,".tirp")  
  file.remove(fname)
}

print("Generating")
curcell <- 0
while(curcell < 10000) {
  rand_sample <- sample(1:nrow(df_align_bc_dist),size = 1)
  for_species <- stringr::str_split_fixed(df_align_bc_dist$species[rand_sample]," ",3)[1]
  
  rand_sample <- sample(1:nrow(df_align_bc_dist),size = 1)
  
  #take from empirical distribution
  num_read <- df_align_bc_dist$cnt[rand_sample]
  
  for_genome <- which(str_starts(names(allfa_main),for_species))
  if(length(for_genome)>0){
    curcell <- curcell + 1
    
    print(paste(curcell,for_species,num_read))
    
    for_genome <- for_genome[1]
    
    out_rand <- round(runif(1, 1,num_shard))
    fname <- paste0("/husky/henriksson/atrandi/simulated2/filtered.",out_rand,".tirp")  
    
    df <- simulate_one_cell(
      as.character(allfa[for_genome]), 
      genome_name = paste0(names(allfa)[for_genome], "#", curcell),
      num_read = num_read,
      insert_max_size = insert_max_size
    )
    write.table(df, fname, row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append=TRUE)   
    
  } else {
    print(paste("no",for_species))
  }
  
}

print("Finalizing")
for(i in 1:num_shard){
  print(i)
  fname <- paste0("/husky/henriksson/atrandi/simulated2/filtered.",i,".tirp")  
  finalize_file(fname)
}
