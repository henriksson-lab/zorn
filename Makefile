build:
	Rscript -e "devtools::document()"
	R CMD build .

build_site:
	Rscript -e "pkgdown::build_site()"

install:
	Rscript -e "install.packages(\"Zorn_0.1.0.tar.gz\", repos = NULL, type = 'source')"
	#install.packages("Zorn_0.1.0.tar.gz", repos = NULL, type = 'source')

addgit:
	git add man
	git add R
	git add vignettes

loc:
	wc -l R/*.R tutorial/*.R tutorial/*.sh


upload_site:
	#

clean:
	rm *.out
