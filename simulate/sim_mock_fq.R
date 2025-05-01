library(Biostrings)


############### Load fasta to simulate
allfa <- Biostrings::readDNAStringSet("/husky/fromsequencer/240809_novaseq_wgs1/trimmed/ref10/all.fa")
data.frame(
  name=names(allfa),
  len=str_length(allfa)
)
allfa_main <- allfa[str_length(allfa)>1000000]



#####################################
##################################### Barcoding and adapters
#####################################

list_bc <- read.csv("atrandi_barcodes.tsv",sep="\t")
rownames(list_bc) <- list_bc$well

generate_r2_prepend <- function() {
  
  bc1 <- sample(list_bc$well[list_bc$pos=="1"],size = 1)
  bc2 <- sample(list_bc$well[list_bc$pos=="2"],size = 1)
  bc3 <- sample(list_bc$well[list_bc$pos=="3"],size = 1)
  bc4 <- sample(list_bc$well[list_bc$pos=="4"],size = 1)
  
  total_bc <- paste0(
    bc1,"_",bc2,"_",bc3,"_",bc4
  )
  list_bc[bc1,]$seq
  
  to_prepend <- paste0(
    list_bc[bc4,]$seq,
    "AGGA",
    list_bc[bc3,]$seq,
    "ACTC",
    list_bc[bc2,]$seq,
    "AAGG",
    list_bc[bc1,]$seq,
    "T"
  )
  
  c(total_bc, to_prepend)
}
#generate_r2_prepend()


#For simplicity, assume that the insert is long enough that we never see the adapters


#####################################
#####################################
#####################################

qread <- str_flatten(rep("F",150))


### TODO: store as FASTQ
#ch = as.character(allfa[1])

simulate_one_cell_fq <- function(ch, genome_name, num_read=1000, insert_max_size=800, to_prepend) {
  totlen <- str_length(ch)
  #insert_max_size <- 800
  start_pos <- round(runif(num_read, min=1, max=totlen-insert_max_size))
  end_pos <- start_pos + round(runif(num_read, min=400, max=insert_max_size))
  
  all_inserts <- str_sub(ch, start_pos, end_pos)
  
  insert_r1 <- str_sub(all_inserts, 1, 150)
  insert_r2_rc <- str_sub(
    all_inserts, 
    str_length(all_inserts)-(149-stringr::str_length(to_prepend)), 
    str_length(all_inserts)
  )
  insert_r2 <- paste0(
    to_prepend, 
    as.character(reverseComplement(DNAStringSet(insert_r2_rc)))
  )

  #Interleave reads for R1
  list_r1 <- c(rbind(
    paste0("@",genome_name,"/1"),
    insert_r1,
    "+",
    qread
  ))
  
  #Interleave reads for R2
  list_r2 <- c(rbind(
    paste0("@",genome_name,"/2"),
    insert_r2,
    "+",
    qread
  ))
  
  list(
    r1=list_r1,
    r2=list_r2
  )  
}


simulate_all_cells_fq <- function(allfa, outf_r1, outf_r2, num_cells=1000, num_read=1000, insert_max_size=800, append=FALSE){
  if(!append){
    file.remove(outf_r1)
    file.remove(outf_r2)
  }
  for(cur_cell in 1:num_cells){
    print(cur_cell)
    for_genome <- runif(1, min = 1, max=length(allfa))
    
    
    to_prepend <- generate_r2_prepend()
    to_prepend_id <- to_prepend[1]
    to_prepend_seq <- to_prepend[2]
    
    rp <- simulate_one_cell_fq(
      as.character(allfa[for_genome]), 
      genome_name = paste0(str_split_i(names(allfa),"\\.",1)[for_genome], "#", to_prepend_id),
      num_read = num_read,
      insert_max_size = insert_max_size,
      to_prepend
    )
      
    writeLines(con = outf_r1, text = rp$r1)
    writeLines(con = outf_r2, text = rp$r2)
  }
}


finalize_file <- function(outf_r1, outf_r2){
  system(paste("pigz -f",outf_r1))
  system(paste("pigz -f",outf_r2))
}


outf_r1 <- "/husky/henriksson/atrandi/simulated3/raw/sim3_R1_001.fastq"
outf_r2 <- "/husky/henriksson/atrandi/simulated3/raw/sim3_R2_001.fastq"

simulate_all_cells_fq(
  allfa,
  outf_r1,
  outf_r2,
  num_read=10000,
  num_cells=100, 
  insert_max_size=800,
  append=FALSE
)

finalize_file(
  outf_r1,
  outf_r2
)




