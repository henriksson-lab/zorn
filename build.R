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
