build:
	Rscript -e "devtools::document()"
	R CMD build .

install:
	Rscript -e "install.packages(\"Zorn_0.1.0.tar.gz\", repos = NULL, type = 'source')"
	#install.packages("Zorn_0.1.0.tar.gz", repos = NULL, type = 'source')


loc:
	wc -l R/*.R tutorial/*.R tutorial/*.sh


upload_site:
	scp 
