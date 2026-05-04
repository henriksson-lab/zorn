build:
	Rscript -e "devtools::document()"
	R CMD build .

build_site:
	Rscript -e "pkgdown::build_site()"

install: build
	Rscript -e "install.packages(\"Zorn_0.5.18.tar.gz\", repos = NULL, type = 'source')"

gitaddall:
	git add man
	git add R/*R
	git add vignettes

loc:
	wc -l R/*.R tutorial/*.R tutorial/*.sh vignettes/*.Rmd


upload_site:
	#

clean:
	rm *.out
