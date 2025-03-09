library(Zorn)

if(FALSE){
  source("R/job_general.R")
  source("R/job_local.R")
  source("R/job_slurm.R")
  source("R/bascet_file.R")
  source("R/zorn.R")
  source("R/shell.R")
  source("R/zorn_aggr.R")
  source("R/count_kmer.R")
  source("R/refgenome.R")
  source("R/kraken.R")
  source("R/container.R")
  
}



bascet_inst <- getBascetImageInstance()

TestBascetInstance(bascet_inst)




#######

bascetRoot <- "/home/m/mahogny/mystore/atrandi/wgs_miseq2/"



#"/home/m/mahogny/mystore/atrandi/wgs_miseq2/con"
#TODO test quast 


rawmeta <- DetectRawFileMeta("/home/m/mahogny/mystore/dataset/atrandi/240903_wgs_atcc2_miseq")

bascet_runner <- LocalRunner(direct = TRUE, show_script=TRUE)


### Debarcode the reads, then sort them.
BascetGetRaw(
  bascetRoot,
  rawmeta,
  runner=bascet_runner,
  bascet_instance = bascet_inst
)




bascetRoot <- "/home/m/mahogny/mystore/atrandi/minitest/"

###
BascetMapCell(
  bascetRoot,
  withfunction = "_quast",
  inputName = "contigs",  ### quast need to run on contigs!!
  outputName = "quast",
  runner=bascet_runner,
  bascet_instance = bascet_inst
)
#--min-contig  can be lowered for quast TODO. make it 100?

### Get QUAST output
quast_aggr <- MapListAsDataFrame(BascetAggregateMap(
  bascetRoot,
  "quast",
  aggr.quast,
  bascet_instance=bascet_instance  
))





if(FALSE){
  
  fquast <- OpenBascet(bascetRoot, "quast")
  allf <- BascetListFilesForCell(fquast, "D1_B4_F7_B12")               ## todo for all cells, and one cell!
  allf
  
  allf[allf$file!="cellmap.log",]
  
  tfile <- BascetReadFile(fskesa, "E2_B4_E9_E11", "transposed_report.tsv", as="tempfile")
  readLines(tfile)
  
  
}

