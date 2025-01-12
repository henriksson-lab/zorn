source("R/job_general.R")
source("R/job_local.R")
source("R/job_slurm.R")
source("R/bascet_file.R")
source("R/zorn.R")
source("R/shell.R")
source("R/zorn_aggr.R")
source("R/count_kmer.R")





inst <- LocalInstance(direct = TRUE, show_script=TRUE)
bascetRoot = "/husky/henriksson/atrandi/wgs_miseq2/"
rawmeta <- DetectRawFileMeta("/husky/fromsequencer/240903_wgs_atcc2_miseq/raw")


### Debarcode the reads, then sort them.
BascetGetRawAtrandiWGS(
  bascetRoot,
  rawmeta,
  runner=inst
)


### Decide cells to include
h <- ReadHistogram(bascetRoot,"debarcoded")
PlotHistogram(h)
includeCells <- h$cellid[h$count>500]   ### for miseq #2
length(includeCells)

### Shardify i.e. divide into multiple sets of files for parallel processing
BascetShardify(
  bascetRoot,
  includeCells = includeCells,
  runner = inst
)

### Assemble ; 
if(FALSE){
  # assembly - 
  #test_skesa:
  #  rm -Rf temp; cargo +nightly run mapcell -i testdata/out_complete.0.tirp.gz -o testdata/skesa.0.zip -s _skesa --show-script-output
  #test_kmc_reads:
}
#  rm -Rf temp; cargo +nightly run mapcell -i testdata/filtered.0.tirp.gz -o testdata/kmc.0.zip -s _kmc_process_reads --show-script-output
#rm -Rf temp; cargo +nightly run mapcell -i testdata/out_complete.0.tirp.gz -o testdata/kmc.0.zip -s _kmc_process_contigs --show-script-output



#"/husky/henriksson/atrandi/wgs_miseq2/"

### Generate KMER databases for each cell
BascetMapCell(
  bascetRoot,
  withfunction = "_kmc_process_reads",
  inputName = "filtered",
  outputName = "kmc",
  runner=inst
)


### Sum up KMC KMER databases
BascetFeaturise(
    bascetRoot, 
    runner=inst
)
getwd()



bascet_instance.default  #temp dir here





