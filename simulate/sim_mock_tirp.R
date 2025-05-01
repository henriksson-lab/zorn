library(Biostrings)


simulate_one_cell <- function(ch, genome_name, num_read=1000, insert_max_size=800) {
  totlen <- str_length(ch)
  #insert_max_size <- 800
  start_pos <- round(runif(num_read, min=1, max=totlen-insert_max_size))
  end_pos <- start_pos + round(runif(num_read, min=400, max=insert_max_size))
  
  all_inserts <- str_sub(ch, start_pos, end_pos)
  
  insert_r1 <- str_sub(all_inserts, 1, 150)
  insert_r2_rc <- str_sub(all_inserts, str_length(all_inserts)-149, str_length(all_inserts))
  insert_r2 <- as.character(reverseComplement(DNAStringSet(insert_r2_rc)))
  
  str_length(insert_r1)
  str_length(insert_r2)
  
  qread <- str_flatten(rep("F",150))
  
  data.frame(
    name=genome_name, #paste0(,"__",1:length(insert_r1)),
    pos1=1,
    pos2=1,
    r1=insert_r1,
    r2=insert_r2,
    q1=qread,
    q2=qread,
    umi=""
  )  
}


simulate_all_cells <- function(allfa, outf, num_cells=1000, num_read=1000, insert_max_size=800, append=FALSE, finalize=TRUE, startseq=1){
  if(!append){
    file.remove(outf)
  }
  for(cur_cell in 1:num_cells){
    print(cur_cell)
    for_genome <- runif(1, min = 1, max=length(allfa))
    
    df <- simulate_one_cell(
      as.character(allfa[for_genome]), 
      genome_name = paste0(str_split_i(names(allfa),"\\.",1)[for_genome], "#", startseq+cur_cell),
      num_read = num_read,
      insert_max_size = insert_max_size
    )
    
    write.table(df, outf, row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append=TRUE)   
  }
  if(finalize){
    print("finalizing")
    finalize_file(outf)
  }
}


finalize_file <- function(outf){
  system(paste0("bgzip -f ",outf))
  system(paste0("tabix -f -p bed ",outf, ".gz"))  
}


############### Load fasta to simulate
allfa <- Biostrings::readDNAStringSet("/husky/fromsequencer/240809_novaseq_wgs1/trimmed/ref10/all.fa")
data.frame(
  name=names(allfa),
  len=str_length(allfa)
)

simulate_all_cells(
  allfa, 
  "/husky/henriksson/atrandi/simulated1/filtered.1.tirp",
  num_cells=1000,
  num_read=1000
)
#1M reads


simulate_all_cells(
  allfa, 
  "/husky/henriksson/atrandi/simulated1/filtered.1.tirp",
  num_cells=100,
  num_read=10000
)
#1M reads total



#####################################
#####################################
#####################################

allfa_main <- allfa[str_length(allfa)>1000000]

simulate_all_cells(
  allfa_main, 
  "/husky/henriksson/atrandi/simulated1/filtered.1.tirp",
  num_cells=30,
  num_read=5000,
  startseq=1000
#  finalize = FALSE
)
simulate_all_cells(
  allfa_main, 
  "/husky/henriksson/atrandi/simulated1/filtered.2.tirp",
  num_cells=30,
  num_read=10000,
  startseq=2000
#  append = TRUE
#  finalize = FALSE
)
simulate_all_cells(
  allfa_main, 
  "/husky/henriksson/atrandi/simulated1/filtered.3.tirp",
  num_cells=30,
  num_read=30000,
  startseq=3000
#  append = TRUE,
#  finalize = FALSE
)
simulate_all_cells(
  allfa_main, 
  "/husky/henriksson/atrandi/simulated1/filtered.4.tirp",
  num_cells=10,
  num_read=50000,
  startseq=4000
#  append = TRUE
)





######################## 
######################## write individual genomes as cells
######################## 

#str_length(allfa_main)

df <- data.frame(
  name=names(allfa_main),
  pos1=1,
  pos2=1,
  r1=as.character(allfa_main),
  r2="",
  q1="",
  q2="",
  umi=""
)
for(i in 1:nrow(df)){
  df$q1 <- str_flatten(rep("F",str_length(df$r1[i])))
}
outf <- "/husky/henriksson/atrandi/simulated1/allgenome/filtered.1.tirp"
write.table(df, outf, row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
finalize_file(outf)



