
if(FALSE){
  source("R/job_general.R")
  source("R/job_local.R")
  source("R/job_slurm.R")
  
  source("R/bascet_file.R")
  source("R/zorn.R")
  source("R/shell.R")
  source("R/zorn_aggr.R")
  source("R/aggr_functions.R")
  source("R/count_kmer.R")
  source("R/countsketch.R")
  source("R/refgenome.R")
  source("R/kraken.R")
  source("R/container.R")
  source("R/ext_tools.R")
  
} else {
  library(Zorn)
}


################################################################################
################ Building documentation website ################################
################################################################################
if(FALSE){
  install.packages("pkgdown")

  # These deps are needed to build the site 
  install.packages("cachem")
  install.packages("fastmap")
  
  # Run once to configure your package to use and deploy pkgdown
  usethis::use_pkgdown_github_pages()  
}


################################################################################
################ Building the package ##########################################
################################################################################

if(FALSE){

  # Preview your site locally before publishing
  pkgdown::build_site()
  
  devtools::document()
  system("R CMD build .")  
  install.packages("Zorn_0.1.0.tar.gz", repos = NULL, type = 'source')

  ### testing
  ??BarnyardPlotMatrix  
  ??OpenBascet

  citation("Zorn") #does not work   .. Error in meta$Priority : $ operator is invalid for atomic vectors

}
