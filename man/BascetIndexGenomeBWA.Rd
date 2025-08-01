% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refgenome.R
\name{BascetIndexGenomeBWA}
\alias{BascetIndexGenomeBWA}
\title{Index a genome using BWA such that it can be used for alignment}
\usage{
BascetIndexGenomeBWA(
  bascetRoot,
  genomeFile,
  overwrite = FALSE,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
}
\arguments{
\item{genomeFile}{Name of FASTA file}
}
\description{
Index a genome using BWA such that it can be used for alignment
}
