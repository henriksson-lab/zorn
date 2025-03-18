################################################################################
################ Building documentation website ################################
################################################################################
if(FALSE){
  install.packages("pkgdown")

  # These deps are needed for build site 
  install.packages("cachem")
  install.packages("fastmap")
}


if(FALSE){
  
  # Run once to configure your package to use and deploy pkgdown
  usethis::use_pkgdown_github_pages()
  
  # Preview your site locally before publishing
  pkgdown::build_site()

  citation("Zorn") #does not work   .. Error in meta$Priority : $ operator is invalid for atomic vectors
}


################################################################################
################ Building the package ##########################################
################################################################################

if(FALSE){
  
  devtools::document()

  system("R CMD build .")
  
  install.packages("Zorn_0.1.0.tar.gz", repos = NULL, type = 'source')
  ??BarnyardPlotMatrix  
  ??OpenBascet
  
  
  
  
}
